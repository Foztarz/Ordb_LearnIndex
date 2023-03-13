# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 03 08
#     MODIFIED:	James Foster              DATE: 2023 03 13
#
#  DESCRIPTION: Loads a set of ordbetareg models and the associated data and
#               functions. The models are then averaged across a range of 
#               plausible model formulas.
#               
#       INPUTS: A folder containing the original data and
#               output from LearnIndex_ordbeta_multimodel.R.
#               
#      OUTPUTS: Best models chosen and the list of imputed models with that formula
#
#	   CHANGES: - save as much as possible
#             - abandon do.call approach (doesn't find object name)
#
#   REFERENCES: Kubinec, R. (2022). 
#               Ordered Beta Regression: A Parsimonious, Well-Fitting Model for 
#               Continuous Data with Lower and Upper Bounds. 
#               Political Analysis, 1-18. doi:10.1017/pan.2022.20
# 
#              Yao Y, Vehtari A, Simpson D & Gelman A (2018). 
#              Using Stacking to Average Bayesian Predictive Distributions (with Discussion). 
#               Bayesian Analysis 13, 917–1007.
# 
#       USAGE:  
#TODO   ---------------------------------------------
#TODO   
#- fix missing chains +
#- chose models of interest +
#- save averaged model

# Starting parameters -----------------------------------------------------
n_iter = 1e3 # number of modelling iterations to run, 1e2 is faster, 1e4 is more accurate
n_pred = 2e1 # number of predictions levels to run, 5e0 is faster, 2e1 is smoother
n_cores = parallel::detectCores()-1 # number of CPUs to use, leave one for other user functions
n_chains =  min(parallel::detectCores()-1, 4) # number of chains to run, 4 is usually sufficient
spare_cpu = TRUE # don't maximise parallelisation, takes longer but preserves performance for other things
ad_delta = 0.95 # closer to 1.0 means higher resolution sampling
use_mice = FALSE # use functions for multiple imputation (now done manually)
# load_loo = TRUE # load calculated LOO for each model previously
save_loo = TRUE #to avoid overwriting, only save if not loading

# . Load packages ---------------------------------------------------------
require(cmdstanr)#most efficient connection to Stan sampler
require(brms)#Bayesian Regression Models using 'Stan', set of tools to build models
#Check 
require(ordbetareg)#load the package for Ordered Beta regression
if(use_mice){require(mice)}#load the package for multiple imputation in case needed


# . Useful functions ------------------------------------------------------

#general purpose
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


# coefficient of variation
CoefVar = function(x){sd(x)/mean(x)} #https://en.wikipedia.org/wiki/Coefficient_of_variation

# . . Specific functions for model comparison -----------------------------

#wrapper for loading and combining brm_multiple model lists
ModComb = function(mlst, ...)
{ brms::combine_models(mlist = mlst, ...) }

#set up functions for looping through datasets

#function performs LOO with imputed dataset _m_
LooM = function(m, model, impdata, ...)
{
  brms::loo(model[[m]],
            newdata = switch(EXPR = class(impdata), 
                             mids = complete(impdata, action = m), # if mice obj, use mice methods 
                             list = impdata[[m]]), # if a list of datasets, extract the dataset
            ...
  )$estimates
}

#function performs LOO for all imputed datasets for each model
Limp = function(model, impdata, ...)
{
  lapply(X = 1:switch(EXPR = class(impdata), 
                      mids = impdata$m, # if mice obj, use mice methods 
                      list = length(impdata)
  ), 
  FUN = LooM,
  model = model,
  impdata = impdata,
  pointwise = TRUE
  )
}

#aggregate list elements of the same dimensions
AggList =  function(lst, fn = median, ...)
{
  ii = rownames(lst[[1]]) # vector of row identifiers
  jj = colnames(lst[[1]]) # vector of column identifiers
  ar = array(data = unlist(lst), 
             dim = c(dim(lst[[1]]), length(lst))) # array with 3rd dimension stacked across list
  ap = apply(X = ar, #  with this array
             MARGIN = 1:2, # apply function across each element in dimensions 1 and 2
             FUN = fn, #function to be applied
             simplify = FALSE,  #output as list
             ...)
  rownames(ap) = ii # same row identifiers as original
  colnames(ap) = jj # same column identifiers as original
  return(ap)
}


# Select and load dataset of interest -------------------------------------
#guess data location in current folder (or working directory when running interactive session)
path_file = tryCatch(expr = file.path(dirname(sys.frame(1)$ofile), "dataR2.csv"),
                     error = function(e)
                     {
                       file.path(getwd(),"dataR2.csv")
                     }
)
if(!file.exists(path_file))
{path_file = FT_select_file(file_type = 'dataR2.csv')}
#TODO check if this is used
# df = read.table(file = path_file,
#                 sep = ',', # ',' for '0.0' decimals csv, ';' for '0,0' 
#                 header = TRUE
# )
## View(df) # hopefully we have loaded the correct dataset
#convert from learning index to proportion correct
# df = within(df, {perc_corr = (LI + 1)/2})
# df = within(df, {Treatment = as.factor(Treatment)})
#find imputed data
if(file.exists(paste0(path_file, '_man_imputed.Rdata')))
{
  imp_load=load(paste0(path_file, '_man_imputed.Rdata'))
  message('...loading:\n',paste0(path_file, '_man_imputed.Rdata'))
}else
{ 
  imp_path = FT_select_file(file_type = '_imputed.Rdata')
  imp_load=load(imp_path)
}
imp_data = get(imp_load)#ensure consistent naming
rm(imp_load)#original named variable no longer needed, remove


# Select and load all models  -------------------------------------------
path_mod = list.files(path = dirname(path_file),
                      pattern = "_allmodels.Rdata",
                      full.names = TRUE
)
if(length(path_mod))
{
  path_mod = path_mod[[1]]
  message('...loading:\n',path_mod)
}else
{
  path_mod = FT_select_file(file_type = "_allmodels.Rdata")
}
nload = load(file = path_mod)#N.B. nload will be a list of names of loaded variables
model_listlist = get(nload)#ensure consistent naming
rm(nload)#remove doubled object

#this is a list of models,
#each element of which contains a list of models with the same formula,
#fitted to different imputed datasets
# length(model_listlist) # number of elements
## [1] 9
# sapply(model_listlist, length) # number of models per formula
## dspeed_age_formula dspeed_age_mu_formula        dspeed_formula     
## 100                   100                   100                   
## dspeed_mu_formula           max_formula        max_mu_formula 
## 100                   100                   100 
## null_formula         treat_formula      treat_mu_formula 
## 100                   100                   100 

# Set up ordbetareg custom family -----------------------------------------
# code taken from:
# https://github.com/saudiwin/ordbetareg_pack/blob/master/R/modeling.R

#N.B. as of 20230215 the posterior_epred_ord_beta_reg function works as expected
#try to find modelling.R without asking the user
#guess local directory
path_R = tryCatch(expr = {file.path(dirname(sys.frame(1)$ofile), 'modeling.R')}, # directory containing this script
                  error = function(e)
                  {#if that fails, try to find the "Documents" folder
                    file.path(getwd(), 'modeling.R')
                  }
)
#guess modeling.R is in that directory
if(file.exists(path_R))
{
  message('...loading:\n',path_R)
}else
{
  path_R = FT_select_file(file_type = "modeling.R")
}
#run and load the custom family definition
source(path_R)#N.B. relies on functions loaded from BRMS

ordbeta_params = .load_ordbetareg(phi_reg = "both")


# Update models missing chains --------------------------------------------


#Some models are missing chains!
#This will be a bit of work to fix...

#find how many chains each
chns = lapply(X = model_listlist,
              FUN = function(i)
              {lapply(X = i, FUN  = brms::nchains)}
)


# . Set up parallel cluster -----------------------------------------------

clt = parallel::makeCluster( if(spare_cpu){floor(n_cores/n_chains)}else{n_cores-1},
                             type = 'PSOCK')
#export variables and functions that will be used by the parallel functions
#this is slow and memory intensive, try to minimise the number of variables loaded
# takes <5 seconds
invisible({parallel::clusterEvalQ(cl = clt, expr = {library("brms")})})
    # invisible({parallel::clusterEvalQ(cl = clt, expr = {library("ordbetareg")})})
ReMissCh = function(mlst,
                    imp_data,
                    prms = ordbeta_params,
                    mchain = max(unlist(
                                        sapply(X = mlst,
                                               FUN = function(x){ sapply(x, brms::nchains) })
                                        )),
                    ...)#passed to brm
{
  for(i in 1:length(mlst))
           {
             if(brms::nchains(mlst[[i]])<mchain)
             {
               md = mlst[[i]]
               mlst[[i]] = brm(formula = brms::bf(formula(md), 
                                           family = prms$family),
                        data = imp_data[[i]],
                        ...
               )
             }
           }
          return(mlst)
}


#run in parallel
#takes less than 2 minutes
system.time({
  # model_listlist = lapply(
    # FUN = ReMissCh,
model_listlist = parallel::parLapply(
    cl = clt, # the parallel cluster, might make things faster
    fun = ReMissCh,
    X = model_listlist,#N.B. this will be overwritten
    mchain = n_chains,#max(unlist(chns)),
    imp_data = df_manimp,
    cores = if(spare_cpu){n_chains}else{1}, # multiple models can be fitted in parallel using parLapply
    chains =  n_chains,
    iter = n_iter,
    init = '0',
    backend = 'cmdstanr',
    prms = ordbeta_params,
    stanvars = ordbeta_params$stanvars,
    control = list( adapt_delta = ad_delta ), # closer to 1.0 means higher resolution sampling
    silent = 2, #don't print lots of iteration information
    refresh = 0
  )
})
## user  system elapsed 
## 1.97    1.33   64.75 

chns_after = lapply(X = model_listlist,
              FUN = function(i)
              {lapply(X = i, FUN  = brms::nchains)}
)
summary(unlist(chns_after))

# . Combine across imputation for each formula ----------------------------

comb_listlist = lapply(X = model_listlist,
                       FUN = ModComb, #wrapper for brms::combine_models
                       check_data = FALSE # by design, these models are not fitted to the same data
)
# class(model_listlist[[1]]) # before each element was a list
## [1] "list"
# sapply(comb_listlist, class) #each element now contains a combined model 
## dspeed_age_formula dspeed_age_mu_formula        dspeed_formula     
## "brmsfit"                   "brmsfit"                   "brmsfit"                   
## dspeed_mu_formula           max_formula        max_mu_formula 
## "brmsfit"                   "brmsfit"                   "brmsfit" 
## null_formula         treat_formula      treat_mu_formula 
## "brmsfit"                   "brmsfit"                   "brmsfit" 

sy = summary(sapply(comb_listlist, ndraws)) # all now 200000
if(any(diff(range(sy))))
  {stop('Models do not have the same length\n', sy)}else
  {print(sy)}


# . Tidy up ---------------------------------------------------------------
rm(model_listlist) # not used beyond this point

# Find best model ---------------------------------------------------------
#could take 15 minutes
path_ic = file.path(dirname(path_mod), paste0(basename(path_mod), 'ic_comb.Rdata'))
if(!file.exists(path_ic))
{
  system.time({
    comb_listlist = parallel::parLapply(cl = clt,
                                  X = comb_listlist,
                                  fun = add_criterion, 
                                  criterion = 'loo'
    )
  })
## user  system elapsed 
## 1.69    1.98  821.22    #memory limits may be slowing this down
}else
{load(path_ic)}

if(save_loo | !file.exists(path_ic))
{
  save(comb_listlist, file = path_ic )
}

#close the parallel cluster when no longer needed
parallel::stopCluster(clt)

#extract loo estimates
l_all = sapply(X= comb_listlist,
                FUN = function(x)
                  {loo(x)$estimates[1]} #loo_elpd (higher is better)
               )

#find best model
best_model_name = names(which.max(l_all))
#sort and assess
l_all = sort(l_all, decreasing = TRUE)
l_diff = l_all - l_all[1]
l_null = l_all - l_all['null_mod']

print(  round(cbind(elpd = l_all, best =  l_diff, null = l_null), 2 )  )
##                       elpd  best  null
## dspeed_noint_mod    -237.87  0.00  2.50  # best model has no interaction of Dspeed & Treatment
## dspeed_noint_mu_mod -238.17 -0.30  2.20
## dspeed_mod          -238.38 -0.52  1.98  # interaction of Dspeed & Treatment within 1.0 of best
## dspeed_mu_mod       -238.47 -0.60  1.90  
## dspeed_age_mu_mod   -239.20 -1.33  1.17
## age_noint_mod       -239.63 -1.77  0.73  # within 1.0 of null model
## treat_mu_mod        -240.11 -2.25  0.26  # within 1.0 of null model
## null_mod            -240.37 -2.50  0.00  
## age_noint_mu_mod    -240.41 -2.55 -0.05  # worse than null model
## treat_mod           -240.49 -2.62 -0.12  # worse than null model
## dspeed_age_mod      -240.58 -2.71 -0.21  # worse than null model
## age_mod             -240.74 -2.88 -0.38  # worse than null model
## age_mu_mod          -241.52 -3.65 -1.15  # worse than null model
## two_interact_mu_mod -241.65 -3.78 -1.28  # worse than null model
## max_mu_mod          -243.04 -5.17 -2.67  # worse than null model
## two_interact_mod    -245.01 -7.14 -4.64  # worse than null model
## max_mod             -247.69 -9.83 -7.33  # worse than null model

#could take >3 minutes
path_lc = file.path(dirname(path_mod), paste0(basename(path_mod), 'loo_comp.Rdata'))
if(!file.exists(path_lc))
{
system.time({
  #do.call version hits recursion limits
  ## Error: C stack usage  36171070 is too close to the limit
  lc =
    with(comb_listlist,#use non-vectorised version to avoid recursion overflow
         #N.B. this will need to be changed for every update to the formula list
    loo_compare(x = dspeed_noint_mod,
                dspeed_noint_mod,
                dspeed_noint_mu_mod,
                dspeed_mod,      
                dspeed_mu_mod,  
                dspeed_age_mu_mod,
                age_noint_mod,      
                treat_mu_mod,
                null_mod,   
                age_noint_mu_mod,
                treat_mod,     
                dspeed_age_mod,
                age_mod,         
                age_mu_mod,
                two_interact_mu_mod,
                max_mu_mod,         
                two_interact_mod,   
                max_mod, 
                criterion = 'loo'   )
    )
})
## user  system elapsed 
## 177.8 0.61    179.4
}
if(save_loo | !file.exists(path_lc))
{
  save(lc, file = path_lc )
}


#having some trouble printing output
#possibly too much in memory
#save and try later


# Model averaging ---------------------------------------------------------


# . Calculate weights -----------------------------------------------------

#N.B. requires all models to have exactly the same numbers of draws
#seems to take a long time, almost 30 min!
path_wt = file.path(dirname(path_mod), paste0(basename(path_mod), 'wt_comb.Rdata'))
if(!file.exists(path_wt))
{
  system.time({
    weights_models = #only way to ensure interpretable names is to name the models explicitly
                     with(comb_listlist,
                       brms::loo_model_weights(x = dspeed_noint_mod,
                                              dspeed_noint_mod,
                                              dspeed_noint_mu_mod,
                                              dspeed_mod,      
                                              dspeed_mu_mod,  
                                              dspeed_age_mu_mod,
                                              age_noint_mod,      
                                              treat_mu_mod,
                                              null_mod,   
                                              age_noint_mu_mod,
                                              treat_mod,     
                                              dspeed_age_mod,
                                              age_mod,         
                                              age_mu_mod,
                                              two_interact_mu_mod,
                                              max_mu_mod,         
                                              two_interact_mod,   
                                              max_mod)
                     )
    })
}
print( cbind(round(sort(weights_models,decreasing = TRUE),3)) )

# user  system elapsed 
# 942.79   96.06 1150.18   
if(save_loo | !file.exists(path_wt))
{
  save(weights_models, file = path_wt )
}

#having some trouble printing output

## Method: stacking
## --- ---
##                    weight
## dspeed_mod          0.258 # highest weighting as we would expect
## dspeed_noint_mu_mod 0.257 # without phi effects, no dspeed-treatment interaction suggested
## age_noint_mod       0.163
## dspeed_age_mu_mod   0.153
## age_noint_mu_mod    0.085
## null_mod            0.084
## dspeed_mu_mod       0.000
## dspeed_noint_mod    0.000
## dspeed_noint_mod    0.000
## age_mod             0.000
## two_interact_mu_mod 0.000
## treat_mu_mod        0.000
## treat_mod           0.000
## age_mu_mod          0.000
## dspeed_age_mod      0.000
## max_mu_mod          0.000
## max_mod             0.000
## two_interact_mod    0.000


# . Average fit -----------------------------------------------------------
path_av = file.path(dirname(path_mod), paste0(basename(path_mod), 'av_model.Rdata'))
if(!file.exists(path_av))
{
system.time({
  average_model = 
    with(comb_listlist,
          brms::pp_average(x = dspeed_noint_mod,
                            dspeed_noint_mod,
                            dspeed_noint_mu_mod,
                            dspeed_mod,      
                            dspeed_mu_mod,  
                            dspeed_age_mu_mod,
                            age_noint_mod,      
                            treat_mu_mod,
                            null_mod,   
                            age_noint_mu_mod,
                            treat_mod,     
                            dspeed_age_mod,
                            age_mod,         
                            age_mu_mod,
                            two_interact_mu_mod,
                            max_mu_mod,         
                            two_interact_mod,   
                            max_mod,
                            method = 'fitted',
                            robust = TRUE
                          )
  )
})
# user  system elapsed 
# 136.37   12.58  151.11 
}

if(save_loo | !file.exists(path_av))
{
  save(average_model, file = path_av )
}


# . Average posterior -----------------------------------------------------
path_po = file.path(dirname(path_mod), paste0(basename(path_mod), 'av_posterior.Rdata'))
if(!file.exists(path_po))
{
system.time({
    average_post = 
      with(comb_listlist, 
           brms::posterior_average(x = dspeed_noint_mod,
                                    dspeed_noint_mod,
                                    dspeed_noint_mu_mod,
                                    dspeed_mod,      
                                    dspeed_mu_mod,  
                                    dspeed_age_mu_mod,
                                    age_noint_mod,      
                                    treat_mu_mod,
                                    null_mod,   
                                    age_noint_mu_mod,
                                    treat_mod,     
                                    dspeed_age_mod,
                                    age_mod,         
                                    age_mu_mod,
                                    two_interact_mu_mod,
                                    max_mu_mod,         
                                    two_interact_mod,   
                                    max_mod,
                                    ndraws = 1e3
                            )
    )
  })
}

if(save_loo | !file.exists(path_post))
{
  save(average_post, file = path_post )
}