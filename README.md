# Ordb_LearnIndex
A set of R scripts for fitting [Ordered Beta Regression](https://doi.org/10.1017/pan.2022.20) models
to learning index data from bumblebees in an automated [Y-maze](https://doi.org/10.3389/fphys.2019.00678), 
to assess the effects of glyphosate treatment on learning.

# Usage
### Set up and install all requirements
This analysis was run in [Rstudio](https://posit.co/download/rstudio-desktop/) ```2022.12.0 Build 353``` 
running [R](https://cran.r-project.org/) ```version 4.2.2```.
Once an ```R``` script is open, run the line where the cursor is currently displayed and advance to the end of that expression by pressing ```Ctrl+Enter```. 
To run the current section (without advancing) press ```Ctrl+Alt+T``` (select "Show document outline" or press ```Ctrl+shift+o``` to view the sections). 
To run a whole script, open that script and _source_ it (```Ctrl+Shift+S``` on Windows).
Once ```R``` and ```Rstudio``` are installed, open [**Install_ordbetareg.R**](Install_ordbetareg.R) and run each section, 
following the instructions on installing [```Rtools```](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html) 
[```cmdStanR```](https://mc-stan.org/cmdstanr/), [```brms```](https://paul-buerkner.github.io/brms/) and [```ordbetareg```](https://cran.r-project.org/web/packages/ordbetareg/vignettes/package_introduction.html).

### Fit potential models
To impute missing data and fit models run [**LearnIndex_ordbeta_manual__impute.R**](LearnIndex_ordbeta_manual__impute.R)
and select the ```.csv``` file containing the original data.
The user also needs to select [```modeling.R```](modeling.R) if not found automatically.  
This saves all models in an ```Rdata``` file ("OrdBmodel_allmodels.Rdata") 
as a nested ```list``` object,  in which each element contains the models fitted for one formula, 
and each sub-element contains a single model fitted to one imputed dataset. 
The ```list``` of imputed datasets used for fitting is also saved as an ```Rdata``` file (name ending with "_man_imputed.Rdata").
*N.B.* This requires multiple processes to be used in parallel and may use a lot of CPU and working memory.
Takes ~4h on a Lenovo ThinkPad with an 11th Gen Intel(R) Core(TM) ```i7-1165G7``` @ ```2.80GHz``` .

### Compare models and select best fit
To compare the different potential model formulas, run [**LearnIndex_ordbeta_multimodel_compare.R**](LearnIndex_ordbeta_multimodel_compare.R) 
and select the locations of the raw data ```.csv``` and output from **LearnIndex_ordbeta_manual__impute.R** (if not in the current working directory).
If specified in the _Starting parameters_ section (```pp_plots = TRUE```), [posterior predictive plots](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssa.12378#rssa12378-fig-0006) 
will be saved for all models using 100 draws from each model (combined across all imputed datasets). 
Model comparison is then performed using [leave-one-out cross-validation](https://mc-stan.org/loo/) 
to identify the set of model parameters with the highest predictive performance (lowest LOO information criterion).
This is first performed as a test using only the first imputed dataset (```lop_test```) to get an impression of how models compare.
This is then performed across all datasets and the LOO results averaged (median) across imputed datasets to summarise typical predictive performance.
*N.B.* This requires multiple processes to be used in parallel and may use a lot of CPU and working memory.
The best fit model is then saved as both the [combined model](https://mc-stan.org/rstan/reference/sflist2stanfit.html) ("best_model.Rdata")
and the ```list``` of models for each imputed dataset ("best_model_list.Rdata") 
in the directory containing the output from **LearnIndex_ordbeta_manual__impute.R**. 
Takes ~30min on a Lenovo ThinkPad with an 11th Gen Intel(R) Core(TM) ```i7-1165G7``` @ ```2.80GHz``` .

### Inspect a model
To inspect the best fit model, run [**LearnIndex_ordbeta_inspect.R**](LearnIndex_ordbeta_inspect.R) and select the raw data ```.csv``` 
and output from **LearnIndex_ordbeta_multimodel_compare.R** (if not in the current working directory).
in the _Starting parameters_ section (```use_best = FALSE```), the ```list``` of all models will be loaded and 
the user asked to write (all or part of) the name of the model to inspect.
Chains from the model of interest for the 1st imputed dataset are displayed 
(*N.B.* the model combined across imputed datasets contains more chains than can be reasonably plotted) 
and calculates the summary parameters for the model coefficients, including [Rhat](https://arxiv.org/abs/1903.08008) convergence diagnostics
for each imputed dataset.

Coefficient estimates are converted back to the original scale and plotted to show 
the [discrete cutoff points](https://cran.r-project.org/web/packages/ordbetareg/vignettes/package_introduction.html),
and estimated interactions of ```Treatment```, ```Dspeed``` and ```Age```, saved in the folder from which the model was loaded. 
These predictions are also saved in a series of ```.csv``` files ("[model_name]..._predictions...[.csv]").
Finally, a set of summary statistics are calculated for the pooled estimates of coefficients across all imputed datasets.
These are saved in a ```.csv``` file ("[model_name]__results.csv").

### File manifest

**- Install_ordbetareg.R**
Setup instructions for fitting ordered beta regression in the [Stan](https://mc-stan.org/) Bayesian statistical modelling language.

**- LearnIndex_ordbeta_manual__impute.R**
Impute missing data, and fit and save ordered beta regression models across all model formulas and imputed datasets.

**- LearnIndex_ordbeta_multimodel_compare.R**
Compare the fitted models using [leave-one-out cross-validation](https://mc-stan.org/loo/) 
and save the model with the best predictive performance.

**- LearnIndex_ordbeta_inspect.R**
Inspect the model fit and save predicted effects and plots of predictions alongside the original data.

**- modeling.R**
A set of functions for fitting ordered beta regression models using ```brms::brm_multiple()``` rather than ```ordbetareg::ordbetareg()```. 
Taken from [```ordbetareg_pack```](https://github.com/saudiwin/ordbetareg_pack/blob/master/R/modeling.R) (author [Robert Kubinec](https://github.com/saudiwin)).

**- Spearman_plot.R**
An example method for plotting the Spearman correlation equivalent of a regression line. 
_N.B._ Results follow the intuitive predictions of the correlation, 
but the simulation shows they are a poor reflection of the process that generated the data 
(at least where multiple data-points have the same value).


