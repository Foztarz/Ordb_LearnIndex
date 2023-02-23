# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 02 01
#     MODIFIED:	James Foster              DATE: 2023 02 23
#
#  DESCRIPTION: Loads a previously fitted ordbetareg model and the associated
#               data and functions for inspection
#               
#       INPUTS: A folder containing the original data, define_ord_betareg.R and
#               output from LearnIndex_ordbeta_impute.R.
#               
#      OUTPUTS: Saves model summary
#
#	   CHANGES: - calculate z-scores from odd ratios
#             - use modified ordbetareg functions
#             - account for NaN rhats (produced by flat priors)
#             - print important summary statistics reported in the paper
#             - allow selection of model from model list
#
#   REFERENCES: Kubinec, R. (2022). 
#               Ordered Beta Regression: A Parsimonious, Well-Fitting Model for 
#               Continuous Data with Lower and Upper Bounds. 
#               Political Analysis, 1-18. doi:10.1017/pan.2022.20
# 
#       USAGE:  
#TODO   ---------------------------------------------
#TODO   
#- Z values +
#- Check compatibility with manual imputation +
#- Check compatibility with ordbetareg v. 0.7.0 (20230212)  +
#- Load from model list list (i.e. not the best model)  +
#- Clean up test statistic calculation (very messy)
#- Account for predictions from mean only models

# Starting parameters -----------------------------------------------------
n_iter = 1e3 # number of modelling iterations to run, 1e2 is faster, 1e4 is more accurate
n_pred = 2e1 # number of predictions levels to run, 5e0 is faster, 2e1 is smoother
use_mice = FALSE # use functions for multiple imputation (now done manually)
use_best = FALSE # Inspect the best model. If FALSE, the user will be asked which model

#load packages
require(cmdstanr)#most efficient connection to Stan sampler
require(brms)#Bayesian Regression Models using 'Stan', set of tools to build models
require(ordbetareg)#load the package for Ordered Beta regression
if(use_mice){require(mice)}#load the package for multiple imputation in case needed

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

#Apply function across all elements of a list of identical arrays
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

# Set up ordbetareg custom family -----------------------------------------
# code taken from:
# https://github.com/saudiwin/ordbetareg_pack/blob/master/R/modeling.R

#TODO check this is necessary
path_R = FT_select_file(file_type = "modeling.R")
#run and load the custom family definition
source(path_R)#N.B. relies on functions loaded from BRMS
ordbeta_params = .load_ordbetareg(phi_reg = 'both') # add the functions to this object

# Select and load dataset of interest -------------------------------------
path_file = FT_select_file(file_type = 'dataR2.csv')

#read in file
df = read.table(file = path_file,
                sep = ',', # ',' for '0.0' decimals csv, ';' for '0,0' 
                header = TRUE
)
#find imputed data
if(file.exists(paste0(path_file, '_man_imputed.Rdata')))
{
  message('...loading:\n',paste0(path_file, '_man_imputed.Rdata'))
  imp_load=load(paste0(path_file, '_man_imputed.Rdata'))
}else
{
  imp_path = FT_select_file(file_type = '_imputed.Rdata')
  imp_load=load(imp_path) 
}
imp_data = get(imp_load)#ensure consistent naming
rm(imp_load)#original named variable no longer needed, remove

# View(df)
#convert from learning index to proportion correct
df = within(df, {perc_corr = (LI + 1)/2})
df = within(df, {Treatment = as.factor(Treatment)})

# Select and load model of interest -------------------------------------
if(use_best)
{
  path_mod = FT_select_file(file_type = "best_model.Rdata")
  #load the model
  nload = load(file = path_mod)#N.B. nload will be a list of names of loaded variables
  #find model list across imputed data
  if(file.exists(sub(pattern = '_model',x =  path_mod, replacement =  '_model_list')))
  {
    message('...also loading:\n',paste0(path_file, '_model_list.Rdata'))
    mload=load(sub(pattern = '_model',x =  path_mod, replacement =  '_model_list'))
  }else
  {
    message('No model list found')
    mload = load(FT_select_file(file_type = '.Rdata'))
  }
  #simplify naming
  comb_model = get(nload)#combined model
  list_model = get(mload)#list of models for each imputed dataset
  rm(nload,mload)#remove old named variables
}else
{#find list of all fitted models
    allmod_pathlist = list.files(path = dirname(path_file),
                                pattern =  '_allmodels.Rdata$',
                                full.names = TRUE)
    if( # if there are any files ending with that name
      any(length( allmod_pathlist  ))
    )
    {
      path_mod = allmod_pathlist[[1]]
      message('...also loading:\n', path_mod)
      lload=load(path_mod)
    }else
    {
      message('No model list found')
      path_mod = FT_select_file(file_type = '.Rdata')
      lload = load(path_mod)
    }
    all_models = get(lload); rm(lload)#rename and tidy up different names
    names_models = names(all_models)#find all model names
    while(!exists(x = 'list_model'))
    {
    message('\nPlease write the name of the model to inspect:\n',
            paste('\t',names_models,'\n')
            )
    mod_name = readline(prompt = 'Model to inspect: ')
    if(mod_name==''){stop('No model selected.')}
    mod_ind = grep(pattern = mod_name, #search for the user input
                   x =  names_models)
    if(any(length(mod_ind)))
    {
      found_name = names_models[mod_ind[1]] #select the first matching name
      message('\n"', found_name, '" selected.')
      list_model = all_models[[ found_name ]] #list of models for each imputed dataset
      comb_model = brms::combine_models(mlist = list_model,
                                        check_data = FALSE) #list of models for each imputed dataset
      rm('all_models')#this list is very large, perhaps better to remove it?
    }else
    {message('\n "',mod_name, '" not found in ', basename(path_mod), '!\n')}
    Sys.sleep(1.0)#goes too fast for the user to see the message on some computers
    }
}
    
# General inspections of fitting with brms -------------------------------
plot(list_model[[1]])# all chains should hover near a single value, looking like "hairy caterpillars"
# N.B. This is just the 1st imputed model, the are length(list_model) total imputed models
# step...
# through...
# plots

#look at parameter summary
sm = summary(comb_model, robust = TRUE)# Rhat values should be as close to 1.000 as possible
if(exists('list_model'))
{
sms = lapply(X = list_model,
             FUN = brms::fixef, # does this do the same
             robust = TRUE)
mdns =  AggList(lst = sms,
          fn = median
          )
          
print(mdns,digits = 4)
  }else
  {
  print(sm)
}
#https://mc-stan.org/rstan/reference/Rhat.html
if(exists('list_model'))
{
  rhats = sapply(X = list_model,
                 FUN = brms::rhat)
  stripchart(x =rhats~rownames(rhats), 
             ylim = c(min(c(0.99, rhats), na.rm = TRUE), # flat prior can produce NaN rhats 
                      max(c(1.1, rhats), na.rm = TRUE) # flat prior can produce NaN rhats
                      ),
             # xlim = 1+c(-1,1)*length(sm$fixed$Rhat)/100,
             ylab = 'Rhat',
             main = "Rhat should be close to 1.0 if chain estimates agree",
             vertical = TRUE,
             method = 'jitter',
             jitter = 0.1,
             pch = 20,
             col = adjustcolor('blue', alpha = 0.5),
             las = 2,
             cex.axis = 0.5
  )
  abline(h = c(1, 1.05),
         col= c(1,2),
         lty = c(1,2))
  print(range(t(rhats), na.rm = TRUE))
  ## 1.000943,  1.001342
}else
{
  with(sm$fixed,
       {
  stripchart(x =Rhat~rownames(sm$fixed), 
             ylim = c(min(c(0.99, sm$fixed$Rhat)), 
                      max(c(1.1, sm$fixed$Rhat))),
             # xlim = 1+c(-1,1)*length(sm$fixed$Rhat)/100,
             ylab = 'Rhat',
             main = "Rhat should be close to 1.0 if chain estimates agree",
             vertical = TRUE,
             method = 'jitter',
             jitter = 0.1,
             pch = 20,
             col = adjustcolor('blue', alpha = 0.5),
             las = 2,
             cex.axis = 0.5
             )
       }
  )
  abline(h = c(1, 1.05),
         col= c(1,2),
         lty = c(1,2))
  print(range(sm$fixed$Rhat, na.rm = TRUE))
}

# estimates are on a logistic scale (0,1) --> (-Inf, +Inf)
if(exists('mdns'))
{
  plogis(mdns['Intercept', 'Estimate'][[1]])
}else{plogis(sm$fixed['Intercept', 'Estimate'])} # mean response rate for Dspeed == 0
# unspecified coefficient types are means of m & c in line equation y = mx + c
if(exists('mdns'))
{
  1/plogis(mdns['phi_Intercept', 'Estimate'][[1]])
}else{1/plogis(sm$fixed['phi_Intercept', 'Estimate'])} # relative variance for Dspeed == 0 (<1 n-shaped distribution, >1 u-shaped)
# phi coefficients are effects on 1/variance
# cutzero and cutone are the bounds at which responses tend to become discrete
print(plogis(sm$spec_pars$Estimate) * 2-1)
##  -0.6393242,  0.3812131

#View IQR for fitted relationships
# plot(brms::conditional_effects(comb_model), points = TRUE)



# Inspect model predictions ----------------------------------------------
all_draws = prepare_predictions(comb_model)
if(exists('list_model'))
{ # do we really get all draws from the combined model?
all_all_draws = lapply(X = list_model,
                       FUN = prepare_predictions)
}
cutzero  =  plogis(all_draws$dpars$cutzero)
cutone  = plogis(all_draws$dpars$cutzero + exp(all_draws$dpars$cutone))

pdf_file = file.path(dirname(path_mod),
                     paste0(
                       if(use_best){basename(path_mod)}else
                                    {found_name},
                       '__results_modelled_cutoffs.pdf')
                     )
pdf(file = pdf_file)
with(df, hist(LI, breaks = 1e3, border = 'darkblue', main  = 'modelled cutoffs'))
#shade the "degenerate" regions
abline(v = 2*c(median(cutzero), median(cutone))-1,
       col = adjustcolor(2, alpha.f = 0.9), lty = 1, lwd = 2)

abline(v = 2*c(sample(cutzero,1e2), sample(cutone,1e2))-1,
       col = adjustcolor(2, alpha.f = 0.01), lty = 1, lwd = 2)
dev.off()
shell.exec.OS(pdf_file)

# Extract conditional effects ---------------------------------------------


#by default 100 predictions per continuous variable (but fewer for their interactions)
system.time(
  {
ordb_cond =brms::conditional_effects(comb_model, 
                                     method = 'posterior_epred', # posterior epred not working
                                     cores =  parallel::detectCores()-1)
  }
)#takes <30 seconds per imputation, prepare to wait approx. 3 min
#using posterior_epred takes <1 min

# Save predictions --------------------------------------------------------

SavePredData = function(ind, nms, fpath)
{
  write.table(x = ordb_cond[[ind]],
              file = paste0(fpath, '_', nms[ind], '.csv'),
              sep = ',',
              row.names = FALSE
              
  ) 
}
#loop through and save for each model parameter
pred_file = file.path(dirname(path_mod), 
                      paste0(
                            if(use_best){basename(path_mod)}else
                            {found_name},
                            '_predictions')
                      )
invisible(
  {
    lapply(X = 1:length(ordb_cond), 
           FUN = SavePredData,
           fpath = pred_file,
           nms = sub(names(ordb_cond), pattern = ':', replacement = '_')
    )
  }
)







# Plot predictions --------------------------------------------------------


# . Dspeed:Treatment ------------------------------------------------------
#extract predictions including all model components (interaction)
if(!is.null(ordb_cond$`Dspeed:Treatment`))
{
pred_data = ordb_cond$`Dspeed:Treatment`
#save as a PDF
pdf_file = file.path(dirname(path_mod), 
                     paste0(
                       if(use_best){basename(path_mod)}else
                                   {found_name},
                       '__results_Dspeed_Treatment.pdf'))
pdf(file = pdf_file)#open a PDF
#plot each treatment
with(subset(df, Treatment == 1), plot(x = Dspeed, y = LI, pch = 20, col = 'darkblue', xlab = 'Dspeed', ylab = 'learning index'))
with(subset(df, Treatment == 2), points(x = Dspeed, y = LI, pch = 20, col = adjustcolor('orange', alpha = 0.5)) )
#add maximum, minimum and no-learning lines
abline(h = c(-1,0,1))

#plot total prediction intervals
with(subset(pred_data, Treatment == 2), 
     {
       polygon(x = c(sort(Dspeed), rev(sort(Dspeed))), 
               y = c(lower__[order(Dspeed)],
                     rev(upper__[order(Dspeed)])
               )*2-1, 
               col = adjustcolor('orange', alpha.f = 50/256),
               border = NA,
               lwd = 0.1
       )
     }
)
with(subset(pred_data, Treatment == 1), 
     {
       polygon(x = c(sort(Dspeed), rev(sort(Dspeed))), 
               y = c(lower__[order(Dspeed)],
                     rev(upper__[order(Dspeed)])
               )*2-1, 
               col = adjustcolor('darkblue', alpha.f = 50/256),
               border = NA,
               lwd = 0.1
       )
     }
)
#plot the median prediction lines
with(subset(pred_data, Treatment == 1), 
     lines(x = sort(Dspeed), y = estimate__[order(Dspeed)]*2-1, 
           col = 'darkblue',
           lwd = 5)
)
with(subset(pred_data, Treatment == 2), 
     lines(x = sort(Dspeed), y = estimate__[order(Dspeed)]*2-1, 
           col = 'orange',
           lwd = 5)
)
#overplot raw data for comparison
with(subset(df, Treatment == 1), 
     points(x = Dspeed, y = LI, 
            pch = 21, lwd = 2,
            col = adjustcolor('darkblue', alpha = 0.5), 
            bg = adjustcolor('white', alpha.f = 0.8) 
     )
)
with(subset(df, Treatment == 2), 
     points(x = Dspeed, y = LI, 
            pch = 21, lwd = 2,
            col = adjustcolor('orange', alpha = 0.5), 
            bg = adjustcolor('white', alpha.f = 0.8) 
     ) 
)
#close the PDF to save
dev.off()
shell.exec.OS(pdf_file)
}
# . Dspeed:Age       ------------------------------------------------------
#extract predictions including all model components (interaction)
if(!is.null(ordb_cond$`Age:Treatment`))
{
pred_data = ordb_cond$`Age:Treatment`
pdf_file = file.path(dirname(path_file), 
                     paste0(
                       if(use_best){basename(path_mod)}else
                                    {found_name}, 
                            '__results_Age_Treatment.pdf')
                     )
#open a PDF file
pdf(file = pdf_file)
#plot the data for each treatment
with(subset(df, Treatment == 1), plot(x = Age, y = LI, pch = 20, col = 'darkblue', xlab = 'Age', ylab = 'learning index'))
with(subset(df, Treatment == 2), points(x = Age, y = LI, pch = 20, col = adjustcolor('orange', alpha = 0.5)) )
#add maximum, minimum and no-learning lines
abline(h = c(-1,0,1))

#plot total prediction intervals
with(subset(pred_data, Treatment == 2), 
     {
       polygon(x = c(sort(Age), rev(sort(Age))), 
               y = c(lower__[order(Age)],
                     rev(upper__[order(Age)])
               )*2-1, 
               col = adjustcolor('orange', alpha.f = 50/256),
               border = NA,
               lwd = 0.1
       )
     }
)
with(subset(pred_data, Treatment == 1), 
     {
       polygon(x = c(sort(Age), rev(sort(Age))), 
               y = c(lower__[order(Age)],
                     rev(upper__[order(Age)])
               )*2-1, 
               col = adjustcolor('darkblue', alpha.f = 50/256),
               border = NA,
               lwd = 0.1
       )
     }
)
#add median lines of prediction
with(subset(pred_data, Treatment == 1), 
     lines(x = sort(Age), y = estimate__[order(Age)]*2-1, 
           col = 'darkblue',
           lwd = 5)
)
with(subset(pred_data, Treatment == 2), 
     lines(x = sort(Age), y = estimate__[order(Age)]*2-1, 
           col = 'orange',
           lwd = 5)
)
#overplot raw data for comparison
with(subset(df, Treatment == 1), 
     points(x = Age, y = LI, 
            pch = 21, lwd = 2,
            col = adjustcolor('darkblue', alpha = 0.5), 
            bg = adjustcolor('white', alpha.f = 0.8) 
     )
)
with(subset(df, Treatment == 2), 
     points(x = Age, y = LI, 
            pch = 21, lwd = 2,
            col = adjustcolor('orange', alpha = 0.5), 
            bg = adjustcolor('white', alpha.f = 0.8) 
     ) 
)
#close PDF to save
dev.off()
shell.exec.OS(pdf_file)
}
# Report contrasts --------------------------------------------------------
#collect coefficients
if(exists('list_model'))
{
  fxf = lapply(X = list_model,
               # FUN = function(x){fixef(x, robust = TRUE)})
               FUN = brms::fixef, 
               robust = TRUE)
  ordb_qnt =  AggList(lst = fxf,
           fn = median
             )
}else{ordb_qnt = fixef(comb_model, robust = TRUE)}
ordb_coef = cbind(all_draws$dpars$mu$fe$b,
                  all_draws$dpars$phi$fe$b)

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
# do something similar with the cutoff bounds
ordb_cutoffs = cbind(
  rbind('degenerate zeroes',
        'degenerate ones'),
  t(
    sapply(X = all_draws$dpars[3:4],
           FUN = quantile,
           probs = c(0.5, 0.025, 0.975),
           
    )
  ),
  
  sapply(X = all_draws$dpars[3:4],
         FUN = empirical.p,
         tails = 1 # cutoffs are one side of 0 by definition 
  ),
  rbind(c('', '', '', '', '', '', '', '', '', ''),
        c('', '', '', '', '', '', '', '', '', ''))
)

colnames(ordb_cutoffs) = colnames(test_statistics)
#collect together
results_to_save = rbind(test_statistics,
                        ordb_cutoffs)

# Save contrasts ----------------------------------------------------------
#save as a csv file
csv_file = file.path(dirname(path_mod),
                     paste0(
                       if(use_best){basename(path_mod)}else
                                   {found_name},
                       '__results.csv'))
write.table(x = results_to_save,
            file = csv_file,
            sep = ',',
            row.names = FALSE
)

#open the results

#should open in default program
shell.exec.OS(csv_file)

