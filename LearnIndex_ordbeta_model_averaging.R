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
#- save averaged model +
#- tidy up +

# Starting parameters -----------------------------------------------------
load_loo = FALSE # load calculated LOO for each model previously (if found)?
save_loo = TRUE #should the results be saved if not loading?
overwrite_loo = TRUE #even if it will overwrite existing files?
#modelling parameters (see LearnIndex_ordbeta_manual__impute)
n_iter = 1e3 # number of modelling iterations to run, 1e2 is faster, 1e4 is more accurate
n_pred = 2e1 # number of predictions levels to run, 5e0 is faster, 2e1 is smoother
n_cores = parallel::detectCores()-1 # number of CPUs to use, leave one for other user functions
n_chains =  min(parallel::detectCores()-1, 4) # number of chains to run, 4 is usually sufficient
spare_cpu = TRUE # don't maximise parallelisation, takes longer but preserves performance for other things
ad_delta = 0.95 # closer to 1.0 means higher resolution sampling
use_mice = FALSE # use functions for multiple imputation (now done manually)

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

#select a file
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


# . . Specific functions for model comparison -----------------------------

#wrapper for loading and combining brm_multiple model lists
ModComb = function(mlst, ...)
{ brms::combine_models(mlist = mlst, ...) }

#TODO check is used
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
df = read.table(file = path_file,
                sep = ',', # ',' for '0.0' decimals csv, ';' for '0,0'
                header = TRUE
)
# View(df) # hopefully we have loaded the correct dataset
# convert from learning index to proportion correct
df = within(df, {perc_corr = (LI + 1)/2})
df = within(df, {Treatment = as.factor(Treatment)})
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
#TODO find out why this happens in non-interactive sessions!
if(any(ceiling(diff(range(chns)))))
{
  message('Chains missing in ', path_mod, ', \nrecalculating...')

# . Set up parallel cluster -----------------------------------------------
  
  clt = parallel::makeCluster( if(spare_cpu){floor(n_cores/n_chains)}else{n_cores-1},
                               type = 'PSOCK')
  #export variables and functions that will be used by the parallel functions
  #this is slow and memory intensive, try to minimise the number of variables loaded
  # takes <5 seconds
  invisible({parallel::clusterEvalQ(cl = clt, expr = {library("brms")})})
  
  ReMissCh = function(mlst, # model list
                      imp_data, # imputed datasets
                      prms = ordbeta_params, # ordered beta parameters
                      mchain = max(unlist(
                                          sapply(X = mlst,
                                                 FUN = function(x){ sapply(x, brms::nchains) })
                                          )), # expected chain number (should be the maximum for all)
                      ...)#passed to brm
  {
    for(i in 1:length(mlst))
             {
               if(brms::nchains(mlst[[i]])<mchain)#if chains are missing
               {
                 md = mlst[[i]] # extract the model
                 mlst[[i]] = brm(formula = brms::bf(formula(md), #refit with the formula
                                             family = prms$family), # and custom family
                          data = imp_data[[i]], # using associated dataset
                          ...
                 )
               }
             }
            return(mlst) # return either the refitted or original model 
  }
  
  
  #run in parallel
  #takes less than 2 minutes
  system.time({
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
}
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
if(any(ceiling(diff(range(sy)))))
  {stop('Models do not have the same length\n', sy)}else
  {print(sy)}


# . Tidy up ---------------------------------------------------------------
rm(model_listlist) # not used beyond this point

# Find best model ---------------------------------------------------------
#could take 15 minutes
path_ic = file.path(dirname(path_mod), paste0(basename(path_mod), '_ic_comb.Rdata'))
if(load_loo & !file.exists(path_ic) | save_loo & overwrite_loo)
{
  if(load_loo){message('\n', path_ic, ' not found, recalculating...')}
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

if(save_loo & overwrite_loo | !file.exists(path_ic))
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

path_el = file.path(dirname(path_mod), paste0(basename(path_mod), '_elpd_tab.csv'))
print(  round(cbind(elpd = l_all, best =  l_diff, null = l_null), 2 )  )

##                       elpd  best  null
## dspeed_noint_mod    -237.84   0.00  2.51  # best model has no interaction of Dspeed & Treatment
## dspeed_noint_mu_mod -238.16  -0.32  2.19
## dspeed_mod          -238.38  -0.54  1.97  # interaction of Dspeed & Treatment within 1.0 of best
## dspeed_mu_mod       -238.45  -0.61  1.90
## dspeed_age_mu_mod   -239.71  -1.87  0.64  # within 1.0 of null model
## treat_mu_mod        -240.12  -2.28  0.23  # within 1.0 of null model
## null_mod            -240.35  -2.51  0.00
## treat_mod           -240.51  -2.67 -0.16  # worse than null model
## dspeed_age_mod      -241.24  -3.40 -0.89  # worse than null model
## age_noint_mu_mod    -241.33  -3.49 -0.98  # worse than null model
## two_interact_mu_mod -242.19  -4.35 -1.84  # worse than null model
## age_mu_mod          -242.31  -4.47 -1.96  # worse than null model
## age_noint_mod       -243.55  -5.71 -3.20  # worse than null model
## max_mu_mod          -243.94  -6.10 -3.60  # worse than null model
## two_interact_mod    -245.73  -7.89 -5.38  # worse than null model
## age_mod             -246.19  -8.35 -5.84  # worse than null model
## max_mod             -248.35 -10.51 -8.00  # worse than null model

#save if specified
if(save_loo & overwrite_loo| !file.exists(path_el))
{
tab_elpd = capture.output(cbind(elpd = l_all, best =  l_diff, null = l_null))
tab_elpd = do.call(what = rbind,
                    args = strsplit(tab_elpd, split = '\\s+')
                    )
write.table(x = tab_elpd,
            file = path_el,
            sep = ',', 
            row.names = FALSE,
            col.names = FALSE
            )
}

#could take >3 minutes
path_lc = file.path(dirname(path_mod), paste0(basename(path_mod), '_loo_comp.Rdata'))
if(load_loo & !file.exists(path_lc) | save_loo & overwrite_loo)
{
  if(load_loo){message('\n', path_lc, ' not found, recalculating...')}
system.time({
  #do.call version hits recursion limits
  ## Error: C stack usage  36171070 is too close to the limit
  lc =
    with(comb_listlist,#use non-vectorised version to avoid recursion overflow
         #N.B. this will need to be changed for every update to the formula list
    loo_compare(x = dspeed_noint_mod,
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
}else
{ load(path_lc) }

#save if asked
if(save_loo & overwrite_loo | !file.exists(path_lc))
{
  save(lc, file = path_lc )
  lc_co = do.call(what = rbind,
                  args = strsplit(x = capture.output(lc), 
                                  split = '\\s+')
                  )
  write.table(x = lc_co,
              file = paste0(path_lc,'.csv'),
              sep = ',', 
              row.names = FALSE,
              col.names = FALSE
  )
}



##                     elpd_diff se_diff
## dspeed_noint_mod      0.0       0.0  
## dspeed_noint_mu_mod  -0.3       2.4  
## dspeed_mod           -0.5       2.3  
## dspeed_mu_mod        -0.6       2.6  
## dspeed_age_mu_mod    -1.9       2.6  
## treat_mu_mod         -2.3       2.8  
## null_mod             -2.5       3.4  
## treat_mod            -2.7       2.3  
## dspeed_age_mod       -3.4       2.3  
## age_noint_mu_mod     -3.5       2.9  
## two_interact_mu_mod  -4.3       3.4  
## age_mu_mod           -4.5       2.9  
## age_noint_mod        -5.7       2.5  
## max_mu_mod           -6.1       3.9  
## two_interact_mod     -7.9       3.5  
## age_mod              -8.3       2.8  
## max_mod             -10.5       3.9 


# Model averaging ---------------------------------------------------------


# . Calculate weights -----------------------------------------------------

#N.B. requires all models to have exactly the same numbers of draws
#seems to take a long time, almost 30 min!
path_wt = file.path(dirname(path_mod), paste0(basename(path_mod), '_wt_comb.Rdata'))
if(load_loo & !file.exists(path_wt) | save_loo & overwrite_loo)
{
  if(load_loo){message('\n', path_wt, ' not found, recalculating...')}
  system.time({
    weights_models = #only way to ensure interpretable names is to name the models explicitly
                     with(comb_listlist,
                       brms::loo_model_weights(x = dspeed_noint_mod,
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
}else
{load(path_wt)}

#display calculated model weights
print( cbind(round(sort(weights_models,decreasing = TRUE),3)) )

# user  system elapsed 
# 1313.14  122.83 1626.75    
if(save_loo & overwrite_loo | load_loo & !file.exists(path_wt))
{
  save(weights_models, file = path_wt )
  wm_co = do.call(what = rbind,
                  args = strsplit(x = capture.output(weights_models), 
                                  split = '\\s+')
  )
  write.table(x = wm_co,
              file = paste0(path_wt,'.csv'),
              sep = ',', 
              row.names = FALSE,
              col.names = FALSE
  )
}

#having some trouble printing output

## Method: stacking
## --- ---
##                    weight
## dspeed_mu_mod       0.238 # highest weighting as we would expect, though without phi effects
## dspeed_mod          0.210 # model with phi effects has 2nd highest weighting 
## null_mod            0.162 # model with only mean
## treat_mu_mod        0.141 # model with only effects of treatment on mean
## dspeed_noint_mu_mod 0.123 # model without interaction of Dspeed and treatment, only mean
## dspeed_noint_mod    0.123 # model without interaction of Dspeed and treatment, effects on mean & variance
## treat_mod           0.002  # model with only effects of treatment on mean & variance
## two_interact_mu_mod 0.001
## age_noint_mu_mod    0.000
## dspeed_age_mu_mod   0.000
## max_mu_mod          0.000
## two_interact_mod    0.000
## dspeed_age_mod      0.000
## age_noint_mod       0.000
## age_mod             0.000
## age_mu_mod          0.000
## max_mod             0.000


# . Average fit -----------------------------------------------------------
#Compute posterior predictive draws averaged across models. 
#seems to take a long time, almost 30 min!
path_av = file.path(dirname(path_mod), paste0(basename(path_mod), '_av_model.Rdata'))
if(load_loo & !file.exists(path_av) | save_loo & overwrite_loo)
{
  if(load_loo){message('\n', path_av, ' not found, recalculating...')}
system.time({
  average_model = 
    with(comb_listlist,
          brms::pp_average(x = dspeed_noint_mod,
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
  ## user  system elapsed 
  ## 1302.28  139.39 1599.67
}else
{load(path_av)}

if(save_loo & overwrite_loo | load_loo & !file.exists(path_av))
{
  save(average_model, file = path_av )
}

# . . Plot predictions ----------------------------------------------------

#combine with first imputed dataset for plotting
#(N.B. the imputed Age values in this dataset are _not_
#the same as the original data, but Dspeed & Treatment are)

pltdata = cbind(df, average_model)

if(save_loo & overwrite_loo | 
   load_loo & !file.exists(paste0(path_av, '_plot.csv'))
)
{
write.table(x = pltdata, 
            file = paste0(path_av, '_plot.csv'),
            sep = ',',
            row.names = FALSE)
}

#plot the raw data
with(pltdata,
     plot(x = Dspeed,
     y = LI,
     ylim = c(-1,1),
     pch = 21,
     col = c(gray(level = 0,alpha =  0.5), 
             adjustcolor('orange3', alpha.f = 0.5)
     )[Treatment],
     bg = gray(level = 0.7,alpha =  0.7),
     lwd = 2
     )
)
abline(h = c(-1,0,1))#reference lines

#plot modelled CI
with( subset(pltdata, Treatment == 1),
  polygon(x = c(sort(Dspeed), sort(Dspeed, decreasing = TRUE)),
       y = c(Q2.5[order(Dspeed)], Q97.5[order(Dspeed,decreasing = TRUE)])*2-1,
       col = gray(level = 0,alpha =  0.5),
       border = NA
       )
)
with( subset(pltdata, Treatment == 2),
  polygon(x = c(sort(Dspeed), sort(Dspeed, decreasing = TRUE)),
       y = c(Q2.5[order(Dspeed)], Q97.5[order(Dspeed,decreasing = TRUE)])*2-1,
       col = adjustcolor('orange', alpha.f = 0.5),
       border = NA
       )
)
#plot modelled central tendency
with( subset(pltdata, Treatment == 1),
 lines(x = sort(Dspeed),
       y = Estimate[order(Dspeed)]*2-1,
       col = 'black',
       lwd = 3
       )
)
with( subset(pltdata, Treatment == 2),
 lines(x = sort(Dspeed),
       y = Estimate[order(Dspeed)]*2-1,
       col = 'orange4',
       lwd = 3
       )
)

# . Average posterior -----------------------------------------------------
#seems to take a long time, almost 10 min!
path_po = file.path(dirname(path_mod), paste0(basename(path_mod), '_av_posterior.Rdata'))
if(load_loo & !file.exists(path_po) | save_loo & overwrite_loo)
{
system.time({
    average_post = 
      with(comb_listlist, 
           brms::posterior_average(x = 
                                    dspeed_mod,      #include only models with the effects from dspeed_mod
                                    dspeed_mu_mod,
                                    dspeed_age_mod,
                                    dspeed_age_mu_mod,
                                    two_interact_mod,   
                                    two_interact_mu_mod,
                                    max_mod,
                                   ndraws = brms::ndraws(dspeed_mod) # just use all
                            )
    )
  })
  ## user  system elapsed 
  ## 474.56   48.39  749.33 
 
}else
{load(path_po)}

#save if requested
if(save_loo & overwrite_loo | load_loo & !file.exists(path_po))
{
  save(average_post, file = path_po )
  write.table(x = cbind(names(average_post),
                        t(apply(X = average_post, 
                                MARGIN = 2, 
                                quantile, 
                                prob = c(0.025, 0.5, 0.975) ))
                        ),
              file = paste0(path_po,'_summary.csv'),
              sep = ',', 
              row.names = FALSE,
              col.names = TRUE
  )
}

#print summary of important coefficients
t(apply(X = average_post,
        MARGIN = 2,
        quantile,
        prob = c(0.025, 0.5, 0.975) ))

##                       2.5%        50%       97.5%
## b_Intercept         -0.11850117  0.1099470  0.33840427
## b_Treatment2        -0.50964643 -0.2157055  0.07293404
## b_Dspeed             0.07419586  0.2678890  0.47250710
## b_Treatment2:Dspeed -0.52102805 -0.2135050  0.08522117
## cutzero             -1.90613025 -1.5208700 -1.14934975
## cutone               0.63693283  0.8032105  0.95831010
## lprior              -5.99291125 -5.6232450  0.00000000

# . . Report contrasts --------------------------------------------------------
#collect coefficients
ordb_qnt = t(apply(X = average_post,
                   MARGIN = 2,
                   quantile,
                   prob = c(0.025, 0.5, 0.975) ))
ordb_qnt = data.frame(Q2.5 = ordb_qnt[ ,'2.5%'],
                          Estimate = ordb_qnt[ ,'50%'],
                          Q97.5 = ordb_qnt[ ,'97.5%']
                          )
row.names(ordb_qnt) = sub(row.names(ordb_qnt), pattern = 'b_', replacement = '')
ordb_coef = average_post

#collect and estimate p values from samples
empirical.p = function(x, tails = 2){ifelse(test = median(x)>0,
                                            yes = mean(x<0)*tails, 
                                            no = mean(x>0)*tails)
}

#collect and estimate standard error from samples
sErr = function(x){sd(x) / sqrt(sum(!is.na(x)))}

#collect and estimate z scores from samples
LOR.z = function(lor, Efun = sErr){median(lor) / Efun(lor)}

#collect and estimate z values from samples
P_z = function(z, tails = 2)
{tails*(1-pnorm(abs(z)))}

#one tailed p-values where CI do not overlap with 0
p_coef1 = apply(X = ordb_coef,
                MARGIN = 2,
                FUN = empirical.p,
                tails = 1
)
#two tailed p-values where CI do not overlap with 0 (more conservative)
p_coef2 = apply(X = ordb_coef,
                MARGIN = 2,
                FUN = empirical.p
)

q_coef = apply(X = ordb_coef,
               MARGIN = 2,
               FUN = quantile,
               probs = c(0.025, 0.975)
)
names(p_coef1) = sub(names(p_coef1), pattern = 'b_', replacement = '')

z_coef = apply(X = ordb_coef,
               MARGIN = 2,
               FUN = LOR.z,
               Efun = sd # for this sample size, sd would make more sense
)

p_gr8r_z = sapply(X = z_coef,
                  FUN = P_z,
                  tails = 2 # more conservative
)
p_one_z = sapply(X = z_coef,
                 FUN = P_z,
                 tails = 1 # reasonable based on other measures, a bit hacked
)

#collect in data frame with original estimates
test_statistics = data.frame(coefficient = names(p_coef1),
                             estimate = unlist(ordb_qnt[names(p_coef1), 'Estimate']),
                             lower_quantile = unlist(ordb_qnt[names(p_coef1),'Q2.5']),
                             upper_quantile = unlist(ordb_qnt[names(p_coef1),'Q97.5']),
                             `prop(H0 draws) 1-tailed` = p_coef1,
                             `prop(H0 draws) 2-tailed` = p_coef2,
                             `one-tailed signif` = ifelse(p_coef1 <0.05, yes = '*', no = ''),
                             `two-tailed signif` = ifelse(p_coef2 <0.05, yes = '*', no = ''),
                             `non-zero CI` = ifelse(
                               apply(X = ordb_qnt[names(p_coef1),c('Q2.5','Q97.5')],
                                     MARGIN = 1,
                                     FUN = function(x){abs(sum(sign(unlist(x))))}), 
                               yes = '*', 
                               no = ''),
                             `Odds-ratio` = ifelse(test = ordb_qnt[names(p_coef1), 'Estimate'] >=0,
                                                   yes = exp(unlist(ordb_qnt[names(p_coef1), 'Estimate'])),
                                                   no = -1/exp(unlist(ordb_qnt[names(p_coef1), 'Estimate']))
                             ),
                             `OR-lower` = ifelse(test = ordb_qnt[names(p_coef1), 'Q2.5'] >=0,
                                                 yes = exp(unlist(ordb_qnt[names(p_coef1), 'Q2.5'])),
                                                 no = -1/exp(unlist(ordb_qnt[names(p_coef1), 'Q2.5']))
                             ),
                             `OR-upper` = ifelse(test = ordb_qnt[names(p_coef1), 'Q97.5'] >=0,
                                                 yes = exp(unlist(ordb_qnt[names(p_coef1), 'Q97.5'])),
                                                 no = -1/exp(unlist(ordb_qnt[names(p_coef1), 'Q97.5']))
                             ),
                             `z-score` = z_coef,
                             `p-greater-z` = p_gr8r_z,
                             `p-one_tailed-z` = p_one_z
)

#collect together
results_to_save = test_statistics

# . . Save contrasts ----------------------------------------------------------
#save as a csv file
csv_file = paste0(path_av,
                 '__statistics.csv')
write.table(x = results_to_save,
            file = csv_file,
            sep = ',',
            row.names = FALSE
)

#open the results

#should open in default program
shell.exec.OS(csv_file)
