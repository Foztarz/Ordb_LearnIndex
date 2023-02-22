# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 24 01
#     MODIFIED:	James Foster              DATE: 2023 02 17
#
#  DESCRIPTION: Loads ordbetareg functions and learning index dataset. Multiply
#               imputes datasets to infer missing ages values. Fits and saves the
#               maximal model, including effects of sensitivity (Dspeed) and age.
#               
#       INPUTS: A folder containing the original data, define_ord_betareg.R and
#               output from LearnIndex_ordbeta_impute.R.
#               
#      OUTPUTS: Saves refitted models
#
#	   CHANGES: - Load updated ordbetareg version that works with LOO
#             - Option to try to run >maximum number of threads in parallel
#             - Model with only effects of Age and Treatment
#
#   REFERENCES: Kubinec, R. (2022). 
#               Ordered Beta Regression: A Parsimonious, Well-Fitting Model for 
#               Continuous Data with Lower and Upper Bounds. 
#               Political Analysis, 1-18. doi:10.1017/pan.2022.20
#               
#               Van Buuren, S. (2018). 
#               Flexible Imputation of Missing Data. Second Edition.
#               Chapman & Hall/CRC. Boca Raton, FL. https://stefvanbuuren.name/fimd/
# 
#               Bürkner P-C. (2018) 
#               Advanced Bayesian Multilevel Modeling with the R Package brms.
#               The R Journal, 10(1), 395–411. doi:10.32614/RJ-2018-017.
# 
#       USAGE:  
#TODO   ---------------------------------------------
#TODO   
#- Save brms models with "file = " argument +
#- Fit null model +
#- Manual imputation  +
#- Model with effects of age and treatment only

# Starting parameters -----------------------------------------------------
n_iter = 1e3 # number of modelling iterations to run, 1e2 is faster, 1e3 is more accurate
n_pred = 2e1 # number of predictions levels to run, 5e0 is faster, 2e1 is smoother
n_mimp = 1e2 # number of imputed datasets to generate, 5e0 is faster, 1e2 is smoother
n_cores = parallel::detectCores()-1 # number of CPUs to use, leave one for other user functions
n_chains = min(c(parallel::detectCores()-1, 4)) # 4 chains is sufficient, combined with imputation we have many
spare_cpu = TRUE # don't maximise parallelisation, takes longer but preserves performance for other things
ad_delta = 0.95 # closer to 1.0 means higher resolution sampling

# Useful functions --------------------------------------------------------

#useful functions
#Open file with default program on any OS
# https://stackoverflow.com/a/35044209/3745353
shell.exec.OS = function(x){
  # replacement for shell.exec (doesn't exist on MAC)
  if (exists("shell.exec",where = "package:base"))
  {return(base::shell.exec(x))}else
  {comm <- paste0('open "',x,'"')
  return(system(comm))}
}

FT_select_file = function(file_type = ".csv",
                          sys_win = Sys.info()[['sysname']] == 'Windows')
{
  ##On University computers, use user profile instead of home directory
  if(sys_win){
    #get rid of all the backslashes
    ltp = gsub('\\\\', '/', Sys.getenv('USERPROFILE'))
  }else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
    ltp = Sys.getenv('HOME')#Easier on Mac
  }
  msg = paste('Please select the',#message to display
              '"', file_type,'"',
              'file')
  here_path = tryCatch(expr = #look in the folder containing this file: sys.frame(1)$ofile
                         {file.path(dirname(sys.frame(1)$ofile))},
                       error = function(e)
                       {#if that fails, try to find the "Documents" folder
                         file.path(ltp,'Documents', 
                                   paste0('*',file_type)
                         )
                       }
  )
  # set path to files
  if(sys_win){#choose.files is only available on Windows
    message('\n\n',msg,'\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    path_file  = choose.files(
      default = here_path,#look where the function is stored
      caption = msg
    )
  }else{
    message('\n\n',msg,'\n\n')
    Sys.sleep(0.5)#goes too fast for the user to see the message on some computers
    path_file = file.choose(new=F)
  }
  #show the user the path they have selected
  if(is.null(path_file))
  {stop('No file selected.')}else
  {print(path_file)}
  return(path_file)
}


#load required packages
require(cmdstanr)#most efficient connection to Stan sampler
require(brms)#Bayesian Regression Models using 'Stan', set of tools to build models
require(ordbetareg)#load the package N.B. for latest version remotes::install_github("saudiwin/ordbetareg_pack")

# Set up ordbetareg custom family -----------------------------------------
# code taken from:
# https://github.com/saudiwin/ordbetareg_pack/blob/master/R/modeling.R

#N.B. as of 20230215 the posterior_epred_ord_beta_reg function works as expected
path_R = FT_select_file(file_type = "modeling.R")
#run and load the custom family definition
source(path_R)#N.B. relies on functions loaded from BRMS

ordbeta_params = .load_ordbetareg(phi_reg = "both")

# Select and load dataset of interest -------------------------------------
#dataset updated to include estimates of age range where not known
path_file = FT_select_file(file_type = 'dataR2.csv')

df = read.table(file = path_file,
                sep = ',', # ',' for '0.0' decimals csv, ';' for '0,0' 
                header = TRUE
)

View(df)
#convert from learning index to proportion correct
df = within(df, {perc_corr = (LI + 1)/2})
df = within(df, {Treatment = as.factor(Treatment)})

# Manually impute missing values ------------------------------------------
# in nearly half of all cases Age is not known
summary(df) # where Age is NA, age range is provided as Age_lower, Age_upper

# N.B. Where Age _is_ known (!is.na(Age)) Age_lower & Age_upper == 0, do not impute!

# . Impute missing values -------------------------------------------------
ManImp = #an expression passed to replicate() that will be executed randomly _n_ times
{
    within(df, #return a dataset with the original structure
           {
             Age = ifelse(test = is.na(Age), # 'yes' action
                          yes = round( # N.B. Age is always a whole number of days: with(df, all(Age[!is.na(Age)] %% 1 == 0) )
                                  runif(n = length(Age), # assume a uniform distribution of missing ages
                                        min = Age_lower,
                                        max = Age_upper
                                      )
                                      ),
                          no = Age
                          )
             } )
}

#seed random number generator for reproducibility
set.seed(seed = 19660621)#founding date of University of Konstanz
# generate _n_ datasets with plausible values of Age
df_manimp = replicate(n = n_mimp,
                      expr = ManImp,
                      simplify = FALSE)#export as a list
summary(df_manimp)
summary(df_manimp[[n_mimp]]) #no NAs in Age

#save imputed dataset
save(df_manimp,
     file = paste0(path_file,'_man_imputed.Rdata')
)


# Choose models to fit ----------------------------------------------------

#model with all possible effects
max_mod = bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment + Dspeed + Age + 
    Treatment:Dspeed + Treatment:Age + Dspeed:Age +
    Treatment:Dspeed:Age, # effects on mean
  phi ~ 0 + Intercept + 
    Treatment + Dspeed + Age + 
    Treatment:Dspeed + Treatment:Age + Dspeed:Age +
    Treatment:Dspeed:Age, # effects on 1/variance
  family = ordbeta_params$family # use the custom family defined above
)

#model with only effects on the mean (not the variance)
max_mu_mod =  bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment + Dspeed + Age + 
    Treatment:Dspeed + Treatment:Age + Dspeed:Age +
    Treatment:Dspeed:Age, # effects on mean
  phi ~ 0 +Intercept, # effects on 1/variance
  family = ordbeta_params$family # use the custom family defined above
)

#model with all 2nd order interactions
two_interact_mod = bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment + Dspeed + Age + 
    Treatment:Dspeed + Treatment:Age + Dspeed:Age, # effects on mean
  phi ~ 0 + Intercept + 
    Treatment + Dspeed + Age + 
    Treatment:Dspeed + Treatment:Age + Dspeed:Age, # effects on 1/variance
  family = ordbeta_params$family # use the custom family defined above
)

#model with all  2nd order interactions on the mean
two_interact_mu_mod = bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment + Dspeed + Age + 
    Treatment:Dspeed + Treatment:Age + Dspeed:Age, # effects on mean
  family = ordbeta_params$family # use the custom family defined above
)

#model without any effects (just an average distribution)
null_mod = bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept, 
  phi ~ 0 + Intercept, 
  family = ordbeta_params$family # use the custom family defined above
)

#model with only effects of Dspeed and Treatment and an independent effect of age
dspeed_age_mod =  bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment + Dspeed + Age +
    Treatment:Dspeed, # effects on mean
  phi ~ 0 + Intercept + 
    Treatment + Dspeed + Age + 
    Treatment:Dspeed, # effects on 1/variance
  family = ordbeta_params$family # use the custom family defined above
)

#model with only mean effects of Dspeed and Treatment and an independent effect of age
dspeed_age_mu_mod =  bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment + Dspeed + Age +
    Treatment:Dspeed, # effects on mean
  family = ordbeta_params$family # use the custom family defined above
)

#model with only effects of Dspeed and Treatment
dspeed_mod =  bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment + Dspeed + 
    Treatment:Dspeed, # effects on mean
  phi ~ 0 + Intercept + 
    Treatment + Dspeed + 
    Treatment:Dspeed, # effects on 1/variance
  family = ordbeta_params$family # use the custom family defined above
)

#model with only effects of Dspeed and Treatment on the mean
dspeed_mu_mod =  bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment + Dspeed + 
    Treatment:Dspeed, # effects on mean
  family = ordbeta_params$family # use the custom family defined above
)

#model with only effects of Dspeed and Treatment
age_mod =  bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment + Age + 
    Treatment:Age, # effects on mean
  phi ~ 0 + Intercept + 
    Treatment + Dspeed + 
    Treatment:Dspeed, # effects on 1/variance
  family = ordbeta_params$family # use the custom family defined above
)

#model with only effects of Dspeed and Treatment on the mean
age_mu_mod =  bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment + Age + 
    Treatment:Age, # effects on mean
  family = ordbeta_params$family # use the custom family defined above
)

#model with only effects of Treatment
treat_mod =  bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment, # effects on mean
  phi ~ 0 + Intercept + 
    Treatment, # effects on 1/variance
  family = ordbeta_params$family # use the custom family defined above
)

#model with only effects of Treatment
treat_mu_mod =  bf( # set up a Bayesian model formula
  perc_corr ~ 0 + Intercept + 
    Treatment, # effects on mean
  family = ordbeta_params$family # use the custom family defined above
)

#collect as a list
form_list = mget(ls(pattern = '_mod')) 
summary(form_list)
## Length Class       Mode
## dspeed_age_mod      5      brmsformula list
## dspeed_age_mu_mod   5      brmsformula list
## dspeed_mod          5      brmsformula list
## dspeed_mu_mod       5      brmsformula list
## max_mod             5      brmsformula list
## max_mu_mod          5      brmsformula list
## null_mod            5      brmsformula list
## treat_mod           5      brmsformula list
## treat_mu_mod        5      brmsformula list
## two_interact_mod    5      brmsformula list
## two_interact_mu_mod 5      brmsformula list

# Fit all models ----------------------------------------------------------

#list function for fitting
OrdBetaImpute = function(formula_i,
                         imp_data,
                         prior_i = NULL,
                         df = df,
                         ... #passed to brm_multiple()
)
{
  if(is.null(prior_i))
  {prior_i = get_prior(formula = formula_i, data = imp_data[[1]])}
  
  tryCatch({
    invisible({
      out_model = brm_multiple(formula = formula_i,
                               data = imp_data,
                               prior = prior_i,
                               ...
      )
    })
  })
}

# . Set up parallel cluster -----------------------------------------------
# spare_cpu = FALSE: run lots of models at once, seems fine for memory (a lot of CPU) 
# N.B. the multithreading uses up a a lot of CPU!
clt = parallel::makeCluster( if(spare_cpu){floor(n_cores/n_chains)}else{n_cores-1},
                            type = 'PSOCK')
#export variables and functions that will be used by the parallel functions
#this is slow and memory intensive, try to minimise the number of variables loaded
    # parallel::clusterExport(cl = clt,
    #                         varlist = c(
    #                                       'form_list',
    #                                     'OrdBetaImpute',
                                        # 'df'#,
    #                                     'df_manimp',
    #                                     'n_chains',
    #                                     'n_iter',
    #                                     'ordbeta_params',
    #                                     'get_prior',
    #                                     'brm_multiple'
    #                         )
    # )
# takes <5 seconds
invisible({parallel::clusterEvalQ(cl = clt, expr = {library("brms")})})

# . Fit maximal model -----------------------------------------------------
#Beware, this may use nearly all CPU (memory load currently well managed)
#this could take 4 h!
system.time({
  model_listlist = parallel::parLapply(
    cl = clt, # the parallel cluster, might make things faster
    X = form_list,
    fun = OrdBetaImpute,
    imp_data = df_manimp,
    cores = if(spare_cpu){n_chains}else{1}, # multiple models can be fitted in parallel using parLapply
    # cores = n_chains, # when parallelising, don't use too many cores
    chains =  n_chains,
    iter = n_iter,
    init = '0',
    backend = 'cmdstanr',
    stanvars = ordbeta_params$stanvars,
    control = list( adapt_delta = ad_delta ), # closer to 1.0 means higher resolution sampling
    silent = 2, #don't print lots of iteration information
    refresh = 0,
    combine = FALSE # return a list of fits
  )
})
# user  system elapsed 
# 8.22    3.61 8390.25
#looking at these numbers there is clearly some efficiency to be gained!

#stop the cluster when no longer using it
parallel::stopCluster(cl = clt)

#save list of fitted models
save(model_listlist,
     file = file.path(dirname(path_file),
                      'OrdBmodel_allmodels.Rdata')
)



