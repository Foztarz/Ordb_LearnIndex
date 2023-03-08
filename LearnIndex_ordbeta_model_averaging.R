# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 03 08
#     MODIFIED:	James Foster              DATE: 2023 03 08
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
#	   CHANGES: - 
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
#- fix missing chains
#- chose models of interest
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
  # ndrw = sapply(X = model_listlist, 
  #               FUN = function(i){sapply(i, ndraws)})
  # sht_md = which(ndrw < max(ndrw))
chns = lapply(X = model_listlist,
              FUN = function(i)
              {lapply(X = i, FUN  = nchains)}
)
# system.time({
# for(ii in sht_md)
# {
#   ch = sapply(X = model_listlist[[ii]], nchains)
#   for(jj in which(ch < max(unlist(chns))) )
#   {
#     md_tmp = model_listlist[[ii]][[jj]]
#     model_listlist[[ii]][[jj]] = brm(formula = formula(md_tmp),
#                                      data = df_manimp[[jj]],
#                                      family = ordbeta_params$family, # use the custom family defined above
#                                      cores = if(spare_cpu){n_chains}else{1}, # multiple models can be fitted in parallel using parLapply
#                                      chains =  n_chains,
#                                      iter = n_iter,
#                                      init = '0',
#                                      backend = 'cmdstanr',
#                                      stanvars = ordbeta_params$stanvars,
#                                      control = list( adapt_delta = ad_delta ), # closer to 1.0 means higher resolution sampling
#                                      silent = 2, #don't print lots of iteration information
#                                      refresh = 0
#                                     )
#   }
# }
# })
## user  system elapsed 
## 22.61    5.70  415.64 

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

summary(sapply(comb_listlist, ndraws)) # all now 200000

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
}
if(save_loo | !file.exists(path_ic))
{
  save(comb_listlist, file = path_ic )
}
  
#close the parallel cluster when no longer needed
parallel::stopCluster(clt)

l_all = sapply(X= comb_listlist,
                FUN = function(x)
                  {loo(x)$estimates[1]}
               )
best_model_name = names(which.max(l_all))

system.time({
  # lc = #maybe don't save output
    do.call(what = loo_compare,
               args = within(comb_listlist,
                             {x = get(best_model_name) # 1st argument always needs to be x
                             rm(list = best_model_name)
                             criterion = 'loo'})
              )
})
## user  system elapsed 
## 103.62    0.54  104.61 

# Model averaging ---------------------------------------------------------

#N.B. requires all models to have exactly the same numbers of draws
#seems to take a long time, almost 30 min!
system.time({
weights_models = do.call(what = brms::loo_model_weights,
                          args = within(comb_listlist,
                                        {
                                        x = get(best_model_name) # 1st argument always needs to be x
                                        rm(list = best_model_name)
                                        # model_names = c(best_model_name,
                                        #                 names(comb_listlist)[
                                        #                   !(names(comb_listlist) %in% best_model_name) ]
                                        #                 )
                                        }
                                        )
                        )
})
# user  system elapsed 
# 1354.74  132.61 1662.41 

## Method: stacking
## ------
##   weight
## dspeed_noint_mod    0.000 
## age_mod             0.000 
## age_mu_mod          0.000 
## age_noint_mod       0.157 
## age_noint_mu_mod    0.112 
## dspeed_age_mod      0.000 
## dspeed_age_mu_mod   0.170 
## dspeed_mod          0.238 #gets highest weighting as we would expect, or are the names scrambled?
## dspeed_mu_mod       0.017 
## dspeed_noint_mu_mod 0.219 
## max_mod             0.000 
## max_mu_mod          0.000 
## null_mod            0.086 
## treat_mod           0.000 
## treat_mu_mod        0.000 
## two_interact_mod    0.000 
## two_interact_mu_mod 0.000 


system.time({
  average_model = do.call(what = brms::pp_average,
                          args = within(ic_list,
                                        {
                                          x = get(best_model_name) # 1st argument always needs to be x
                                          rm(list = best_model_name)
                                          method = 'fitted'
                                        }
                          )
  )
})
# user  system elapsed 
# 136.37   12.58  151.11 


system.time({
  average_post = do.call(what = brms::posterior_average,
                          args = within(ic_list,
                                        {
                                          x = get(best_model_name) # 1st argument always needs to be x
                                          rm(list = best_model_name)
                                          ndraws = 1e3
                                        }
                          )
  )
})