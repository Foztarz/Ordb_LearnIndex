# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 02 10
#     MODIFIED:	James Foster              DATE: 2023 02 15
#
#  DESCRIPTION: Loads a set of ordbetareg models and the associated data and
#               functions. The models are then compared using leave-one-out
#               cross validation. Estimates are averaged across imputed datasets.
#               
#       INPUTS: A folder containing the original data, modeling_JJF.R and
#               output from LearnIndex_ordbeta_multimodel.R.
#               
#      OUTPUTS: Best models chosen and the list of imputed models with that formula
#
#	   CHANGES: - save multimodel loo
#
#   REFERENCES: Kubinec, R. (2022). 
#               Ordered Beta Regression: A Parsimonious, Well-Fitting Model for 
#               Continuous Data with Lower and Upper Bounds. 
#               Political Analysis, 1-18. doi:10.1017/pan.2022.20
# 
#               Vehtari, A., Gelman, A., & Gabry J. (2016). 
#               Practical Bayesian model evaluation using leave-one-out 
#               cross-validation and WAIC. 
#               Statistics and Computing, doi:10.1007/s11222-016-9696-4.
# 
#       USAGE:  
#TODO   ---------------------------------------------
#TODO   
#- Load model list  +
#- Combine all models +
#- Model comparison across any model list +
#- Reduce memory usage for parallel processing  +
#- Exclude unused modelling functions +
#- Troubleshoot pp_check_ordbeta

# Starting parameters -----------------------------------------------------
n_iter = 1e3 # number of modelling iterations to run, 1e2 is faster, 1e4 is more accurate
n_pred = 2e1 # number of predictions levels to run, 5e0 is faster, 2e1 is smoother
n_cores = parallel::detectCores()-1 # number of CPUs to use, leave one for other user functions
n_chains =  min(parallel::detectCores()-1, 4) # number of chains to run, 4 is usually sufficient
pp_plots = FALSE # plot and show prior predictive distribution for all models
use_mice = FALSE # use functions for multiple imputation (now done manually)
save_loo = TRUE # save calculated LOO for each model (may overwrite previous saved LOO)

# . Load packages ---------------------------------------------------------
require(cmdstanr)#most efficient connection to Stan sampler
require(brms)#Bayesian Regression Models using 'Stan', set of tools to build models
#Check 
require(ordbetareg)#load the package for Ordered Beta regression
if(use_mice){require(mice)}#load the package for multiple imputation in case needed
if(pp_plots){require(ggplot2)}#load the package for %>% operations if using pp_check


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
df = read.table(file = path_file,
                sep = ',', # ',' for '0.0' decimals csv, ';' for '0,0' 
                header = TRUE
)
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

## View(df) # hopefully we have loaded the correct dataset
#convert from learning index to proportion correct
df = within(df, {perc_corr = (LI + 1)/2})
df = within(df, {Treatment = as.factor(Treatment)})

# Select and load all models  -------------------------------------------
path_mod = FT_select_file(file_type = "_allmodels.Rdata")

nload = load(file = path_mod)#N.B. nload will be a list of names of loaded variables

model_listlist = get(nload)#ensure consistent naming
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



# Posterior predictive check ----------------------------------------------

if(pp_plots) # In CRAN version make sure to set 'reverse_bounds' to FALSE
{
  #open a file to save in
  pdf_file = file.path(dirname(path_file), paste0(nload, '_pp_checks.pdf'))
  pdf(file = pdf_file)#open a PDF file
  #plot subsample of predictions for each model
  for(i in 1:length(comb_listlist))
  {
    #use the ordbetareg function for plotting posterior predictive checks
    plt_tmp = ordbetareg::pp_check_ordbeta(comb_listlist[[i]], 
                                       ndraws = 100, reverse_bounds = FALSE)
    #plot the proportions of discrete and condituous data with the model details
    plot(plt_tmp$discrete + 
      ggtitle(label = names(comb_listlist)[i], # model name as title
      subtitle = formula(comb_listlist[[i]])) # model formula as subtitle
    )
    #plot the contiuous component of the data and predictions
    plot(plt_tmp$continuous+ 
      ggtitle(label = names(comb_listlist)[i],
              subtitle = formula(comb_listlist[[i]]))
    )
  }
  dev.off() # close the PDF file to save
  shell.exec.OS(pdf_file) # view the saved PDF file
}


# Model comparison --------------------------------------------------------


# . Set up parallel cluster -----------------------------------------------
#this will require a lot of processing that can be sped up using multiple CPUs
clt = parallel::makeCluster(n_cores-2, # the multithreading uses up a a lot of CPU!
                            type = 'PSOCK')
#export Bayesian modelling functions to the parallel cluster
# takes <5 seconds
invisible({parallel::clusterEvalQ(cl = clt, expr = {library("brms")})})

# . Test pointwise LOO with one dataset -----------------------------------

#single (imputed) dataset estimates give a good approximation of the full combined model
#takes almost 10 minutes
system.time(
  {
    lop_test = parallel::parLapply(cl = clt,
                         X = comb_listlist,
                      fun = brms::loo,
                     newdata =  switch(EXPR = class(imp_data), 
                                       mids = complete(imp_data, action = 1), 
                                       list = imp_data[[1]]), # first imputed dataset
                     pointwise = TRUE
                    )
  }
)
## user  system elapsed 
## 0.97    0.86  409.52 #mostly overhead from importing and exporting to cluster
if(save_loo)
{
  #TODO allow loading of previously generated test
  save(lop_test,file = dirname(path_file), paste0(nload, 'LOO_test.Rdata'))
}

print(loo_compare(x = lop_test))
##                         elpd_diff se_diff
## dspeed_mod           0.0       0.0   
## dspeed_mu_mod        0.0       3.0   
## dspeed_age_mu_mod   -1.2       3.1   
## dspeed_age_mod      -1.3       1.6   
## treat_mu_mod        -1.7       4.1   
## null_mod            -1.9       4.6   
## treat_mod           -2.0       3.7   
## age_mod             -3.2       3.3   
## two_interact_mu_mod -3.5       3.6   
## age_mu_mod          -3.9       4.3   
## max_mu_mod          -4.6       4.0   
## two_interact_mod    -6.1       2.3   
## max_mod             -8.1       3.6   


# . Run LOO across all imputed datasets -----------------------------------

#set up functions for looping through datasets
#export new functions that will be used by the parallel functions
# should be fast
## user  system elapsed 
## 0.00    0.00    0.36
parallel::clusterExport(cl = clt,
                        varlist = c(
                          'Limp',
                          'LooM'
                        )
)
#include functions for multiple imputation 
if(use_mice)
{ invisible({parallel::clusterEvalQ(cl = clt, expr = {library("mice")})}) }

#apply Leave-One-Out validation to all models for all imputed datasets
#could take up to 15 minutes for 7 models
system.time({
    loom_all = parallel::parLapply(cl = clt,
                         X = model_listlist,
                         fun = Limp,
                         impdata = imp_data
    )
})
# user  system elapsed 
# 1.04    0.78  646.40 
if(save_loo)
{
  #TODO allow loading of previously generated list of LOO model validations
  save(loom_all,file = dirname(path_file), paste0(nload, 'LOO_m_all.Rdata'))
}
#for quick comparison, take the median values across imputed datasets
#should be very fast
loom_agg = parallel::parLapply(cl = clt,
                     X = loom_all,
                     fun = AggList,
                     fn = median
                    )
#check the coefficient of variation across imputed models, this should be _very_ small
loom_cov = parallel::parLapply(cl = clt,
                     X = loom_all,
                     fun = AggList,
                     fn = CoefVar # sd/mean
                    )

#close the parallel cluster when no longer needed
parallel::stopCluster(clt)

#each element of loom_agg is now a LOO summary table with median values
# loom_agg[1]
    ## $dspeed_age_formula
    ## Estimate  SE      
    ## elpd_loo -240.1121 7.60479 
    ## p_loo    12.49652  2.541966
    ## looic    480.2243  15.20958

# imputed datasets are very similar, and so model predictive power should vary little
summary(abs(unlist(loom_cov))) #expect SD <10% of mean
    ## Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## 0.0004901 0.0010619 0.0057453 0.0148961 0.0221784 0.0622748   

# Save the model with lowest LOO information criterion --------------------
#sort LOO criteria lowest (least information loss) to highest (worst information loss)
# loo_ic is the 3rd element (there must be a more elegant way)
Extract_looic_est = function(x){x[[3]]} 
#before using LOO IC for model selection, check variation across imputed datasets was low
summary( sapply(X = loom_cov, FUN = Extract_looic_est) )
  ## Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  ## 0.0005109 0.0006421 0.0008664 0.0009559 0.0009834 0.0020706 

#extract LOO IC for all models
all_looic = sort( 
                  sapply(X = loom_agg,
                   FUN = Extract_looic_est
                   )
                  )
print(t(t(all_looic)))
    ## [,1]
    ## dspeed_mod          476.5723
    ## dspeed_mu_mod       476.9219
    ## dspeed_age_mod      478.9725
    ## dspeed_age_mu_mod   479.3268
    ## treat_mu_mod        480.1619
    ## null_mod            480.7124
    ## treat_mod           480.8620
    ## age_mod             483.1547
    ## two_interact_mu_mod 483.7303
    ## age_mu_mod          484.7352
    ## max_mu_mod          485.8961
    ## two_interact_mod    487.5536
    ## max_mod             491.0355

#ELPD
    ## dspeed_mod          -238.2861
    ## dspeed_mu_mod       -238.4609
    ## dspeed_age_mod      -239.4863
    ## dspeed_age_mu_mod   -239.6634
    ## treat_mu_mod        -240.0810
    ## null_mod            -240.3562
    ## treat_mod           -240.4310
    ## age_mod             -241.5774
    ## two_interact_mu_mod -241.8651
    ## age_mu_mod          -242.3676
    ## max_mu_mod          -242.9480
    ## two_interact_mod    -243.7768
    ## max_mod             -245.5178

#select the model with the lowest LOO information criterion
best_model_name = names(all_looic)[[1]] # this will be first in the sorted list
#find the combined model and the list it was combined from
best_model = comb_listlist[[best_model_name]]
best_model_list = model_listlist[[best_model_name]]
#inspect this model
summary(best_model)
    ## Population-Level Effects: 
    ##   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept                 0.11      0.11    -0.11     0.33 1.00   170308   151175
    ## Treatment2               -0.21      0.15    -0.49     0.08 1.00   183309   155883
    ## Dspeed                    0.27      0.10     0.08     0.47 1.00   183486   154427
    ## Treatment2:Dspeed        -0.24      0.15    -0.54     0.06 1.00   180627   152465
    ## phi_Intercept             1.31      0.17     0.97     1.64 1.00   169948   148710
    ## phi_Treatment2           -0.16      0.25    -0.65     0.32 1.00   162594   153908
    ## phi_Dspeed               -0.10      0.15    -0.40     0.17 1.00   183094   149112
    ## phi_Treatment2:Dspeed    -0.29      0.25    -0.77     0.19 1.00   177658   155249

#save each of these with a recognisible name
save(best_model,
     file = file.path(dirname(path_mod), 'best_model.RData') )
save(best_model_list,
     file = file.path(dirname(path_mod), 'best_model_list.RData') )
                 