# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2022 12 19
#     MODIFIED:	James Foster              DATE: 2023 02 15
#
#  DESCRIPTION: Install and set up packages required for Ordered Beta Regression
#               
#       INPUTS: 
#               
#      OUTPUTS: Test result
#
#	   CHANGES: - install from CRAN rather than GitHub
#
#   REFERENCES: Kubinec, R. (2022). 
#               Ordered Beta Regression: A Parsimonious, Well-Fitting Model for 
#               Continuous Data with Lower and Upper Bounds. 
#               Political Analysis, 1-18. doi:10.1017/pan.2022.20
#               
#               Gabry J, Češnovar R, Johnson A (2022). 
#               cmdstanr: R Interface to 'CmdStan'.
#               https://mc-stan.org/cmdstanr/
# 
#               Bürkner, P.-C. (2018). 
#               Advanced Bayesian Multilevel Modeling with the R Package brms. 
#               The R Journal 10, 395–411.
# 
#               Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., 
#               Betancourt, M., Brubaker, M., Guo, J., Li, P. and Riddell, A. (2017). 
#               Stan: A Probabilistic Programming Language. 
#               Journal of Statistical Software 76 doi: 10.18637/jss.v076.i01
# 
#       USAGE:  
#TODO   ---------------------------------------------
#TODO   

#Open file with default program on any OS
# https://stackoverflow.com/a/35044209/3745353
shell.exec.OS = function(x){
  # replacement for shell.exec (doesn't exist on MAC)
  if (exists("shell.exec",where = "package:base"))
  {return(base::shell.exec(x))}else
  {comm <- paste0('open "',x,'"')
  return(system(comm))}
}

# Install required packages from source -----------------------------------
#first download and install Rtools4.3 (N.B. 4.2 had some bugs!)
shell.exec.OS('https://cran.r-project.org/bin/windows/Rtools/')

#install the package for Bayesian modelling on Windows
install.packages('remotes')#install the package for installing remote packages
remotes::install_github("stan-dev/cmdstanr") #follow all instructions

#install accompanying executable program
require(cmdstanr)
install_cmdstan(cores = parallel::detectCores()-1)# may take a while & print a lot!
#if this fails, try https://github.com/stan-dev/cmdstan/archive/refs/tags/v2.30.1.zip

#install package for writing Bayesian models
install.packages('brms')

#install package for ordered beta modelling
#As of 20230215 it should work from CRAN
# https://github.com/saudiwin/ordbetareg_pack/issues/11#issuecomment-1430837740
install.packages('ordbetareg')

#ordbetareg should be installed from GitHub for latest version
    # remotes::install_github("saudiwin/ordbetareg_pack",
    #                         # build_vignettes=TRUE,
    #                         dependencies = TRUE)
#for more information see https://cran.r-project.org/web/packages/ordbetareg/vignettes/package_introduction.html

#if performing multiple imputation algorithmically install the relevant package
# install_cmdstan('mice')#20230215 the imputation is now manual, no longer used

# Test installed packages -------------------------------------------------

#test cmdstanr
cmdstanr::check_cmdstan_toolchain()
cmdstanr::set_cmdstan_path(path = cmdstanr::cmdstan_path())

N <- 200
Y <- rnorm(N, 5, 2)
test_mod = 
  '
data {
  int<lower = 0> N;
  vector[N] Y;
}
parameters {
  real mu;
  real<lower = 0> sigma;
}
model {
  Y ~ normal(mu, sigma);
}
'
write.table(test_mod, 
            file = file.path(getwd(), 'model1.stan'), 
            col.names = F, 
            row.names = F, 
            quote = F)
# set_cmdstan_path(cmdstanpath)
model1 <- cmdstanr::cmdstan_model("model1.stan", quiet = FALSE)
fit1 <- model1$sample(data = list(N = N, Y = Y),
                      seed = 1, chains = 4,
                      iter_sampling = 1000, iter_warmup = 1000,
                      thin = 1)
fit1$summary()


options(mc.cores = parallel::detectCores()-1)


require(brms)
require(ordbetareg)
require(dplyr)

#test ordbetareg
#example from vigentte
data("pew")

model_data <- select(pew[sample(1:(dim(pew)[1]), size = 2e1, replace = FALSE),
                         ],#small subsample just for testing
                         therm,age="F_AGECAT_FINAL",
                     sex="F_SEX_FINAL",
                     income="F_INCOME_FINAL",
                     ideology="F_IDEO_FINAL",
                     race="F_RACETHN_RECRUITMENT",
                     education="F_EDUCCAT2_FINAL",
                     region="F_CREGION_FINAL",
                     approval="POL1DT_W28",
                     born_again="F_BORN_FINAL",
                     relig="F_RELIG_FINAL",
                     news="NEWS_PLATFORMA_W28") %>% 
  mutate_at(c("race","ideology","income","approval","sex","education","born_again","relig"), function(c) {
    factor(c, exclude=levels(c)[length(levels(c))])
  }) %>% 
  # need to make these ordered factors for BRMS
  mutate(education=ordered(education),
         income=ordered(income))

#run and time
  system.time(
  {
ord_fit_mean <- ordbetareg(formula=
                             bf(therm ~ mo(education)*mo(income) +
                             (1|region),
                             phi ~ mo(income)
                             ),
                           data=model_data,
                           phi_reg = TRUE,
                           cores=2,
                           chains=2,
                           iter=1e3,
                           refresh=0,
                           true_bounds = c(0,100),
                           backend = 'cmdstanr'
                           )
  }
)
#takes <40s (mostly compile time)

summary(ord_fit_mean)#Rhat values close to 1.00

conditional_effects(ord_fit_mean)

# Test validation functions -----------------------------------------------
require(ggplot2)
ordbetareg::pp_check_ordbeta(ord_fit_mean, ndraws = 100)

pe = posterior_epred(ord_fit_mean)
pp = posterior_predict(ord_fit_mean)
ep = sum(log_lik(ord_fit_mean))

#takes <2seconds
waic(ord_fit_mean)

#takes <20 seconds
loo(ord_fit_mean, pointwise = TRUE)

system.time({
add_criterion(ord_fit_mean, c("waic", "loo_R2"))
})

#not used by us
# kfold(ord_fit_mean) # takes <120 seconds
##Error: New factor levels are not allowed. #probably specific to this dataset

#not working
#loo_subsample(ord_fit_mean)
##Error: 'observations' is larger than the total sample size in 'data'.

