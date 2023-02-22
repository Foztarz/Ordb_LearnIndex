# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2023 02 15
#     MODIFIED:	James Foster              DATE: 2023 02 22
#
#  DESCRIPTION: Simulates an example dataset for Spearman correlation, and plots
#               the equivalent regression line for that dataset and a user specified one.
#               
#       INPUTS: Starting parameters for simulation. 
#               A file containing input data as a .csv .
#               
#      OUTPUTS: 
#
#	   CHANGES: - commented
#
#   REFERENCES: Spearman C. (1904).
#               The proof and measurement of association between two things. 
#               American Journal of Psychology. 15 (1): 72–101.
#               doi:10.2307/1412159. JSTOR 1412159.
#               
#       USAGE:  
#TODO   ---------------------------------------------
#TODO   
#- Different methods for transforming from rank(y) -> y


# Simulated example -------------------------------------------------------

set.seed(101) # seed the random number generator for repeatability

# . Simulation parameters -------------------------------------------------

nn = 1e2 # number of observations
cf = 1.0 # scaling coefficient
rng = 4 # range of simulated data
snr = 0.8 # heuristic of signal to noise ratio
err = 0.25 # normally distributed error to add


# . Generate simulated data -----------------------------------------------

#data that is not normally distributed along both axes
nonnorm_data = data.frame(x = runif(nn, # generate uniform distribution that increases linearly
                                    min = seq(from = -rng/2, to  = rng/2, length.out = nn)-rng/snr,
                                    max = seq(from = -rng/2, to  = rng/2, length.out = nn)+rng/snr
                                    ),
                          y = (pnorm( # transform to bounded (normal CDF: "inverse probit")
                                  cf*runif(nn, # original scale is also a uniform distribution that increases linearly
                                          min = seq(from = -rng/2, to  = rng/2, length.out = nn)-rng/snr,
                                          max = seq(from = -rng/2, to  = rng/2, length.out = nn)+rng/snr
                                          ) + rng*err*rnorm(1e2) # add noise
                                        )
                              ),
                          x_mu = seq(from = -rng/2, to  = rng/2, length.out = nn), # population mean x
                          y_mu = pnorm(cf*seq(from = -rng/2, to  = rng/2, length.out = nn)) # population mean y
                          )


# . Inspect simulated data ------------------------------------------------

#a plot shows that the y variable is bounded
with(nonnorm_data, plot(x = x, y = y));abline(h = c(0,1))

#Pearson correlation 'fails' to find significant relationship
print(with(nonnorm_data, cor.test(x, y, method = 'pearson')))
  ##r = 0.1802591  ; p-value = 0.07271

#Spearman rank correlation finds the correlation
print(with(nonnorm_data, cor.test(x, y, method = 'spearman')))
  ##r = 0.2065287 ; p-value = 0.03941

#N.B. Spearman correlation is equivalent to
nonnorm_data = within(nonnorm_data, 
                      {
                      x_rank = rank(x)
                      y_rank = rank(y)
                      }
                      )
with(nonnorm_data, cor(x_rank, y_rank))#correlation on the ranks of each variable
  ## [1] 0.2065287


# . Generalise the linear relationship between ranks ----------------------

#Construct empirical cumulative density functions for x and y
Ecdf_x = with(nonnorm_data, ecdf(x)) # function can be applied to values of x
Ecdf_y = with(nonnorm_data, ecdf(y)) # function can be applied to values of y 

# convert x and y to their expected cumulative density (proportion of smaller values)
nonnorm_data = within(nonnorm_data, 
                      {
                        x_cdf = Ecdf_x(x) # find predicted position in range of x
                        y_cdf = Ecdf_y(y) # find predicted position in range of x
                      }
)
#N.B. empirical cumulative density is identical to rank, except that it is scaled
#between 0 and 1, where 0 is the minimum observed value and 1 is the maximum

#fit a linear model to predicted position (cumulative position in rank between 0 and 1)
lm_rxy = with(nonnorm_data, lm(y_cdf ~ x_cdf))
#here we can see the linear relationship between ranks, where one exists (N.B. too many identical values mask this relationship)
with(nonnorm_data, plot(x = x_cdf, y = y_cdf)) # plot the ECDF
abline(lm_rxy, col = 'red') # add the regression line


# . Predict the relationship implied by the rank correlation --------------

#generate a more continuous x sequence to use for predictions
xx = with(nonnorm_data, seq(from = min(x), to  = max(x), length.out = 1e3))

#predict cumulative y value
pred_rxy = with(nonnorm_data,
                {
                predict(object = lm_rxy, #generate regression predictions
                        newdata = data.frame(x_cdf = Ecdf_x(xx)) #using the predicted cumulative density of x
                        )
                }
                )
#plot the predicted CD of y for the continuous x sequence
plot(xx, pred_rxy) # expected y in range of all y values (0,1)

#predicted y value based on minimum and range of y
pred_y = with(nonnorm_data,
              {
              diff(range(y))* # scale by range (difference between minimum and maximum)
                  (pred_rxy-0.5) + # subtract expected median (CD of 0.5)
                  median(y) # add observed median
              }
                )

#smooth prediction with a spline
spl_y = predict(object = smooth.spline(x = xx, 
                                      y = pred_y, 
                                      nknots = round(sqrt( # round to whole number
                                                length(nonnorm_data$x))
                                                )) # use a small number of knots
                )$y # keep only y prediction


# . Plot the prediction based on ranking (CDF) ----------------------------
with(nonnorm_data, plot(x = x, y = y)) # plot original data
with(nonnorm_data,  # plot the minimum, maximum and median of the data 
    {
       abline(h = c(range(y, na.rm = TRUE),0.5), # (potential tether points for CDF -> y conversion)
            lty = c(1,1,2)) 
     }
     ) # include a dashed line of no correlation
lines(xx, spl_y, col = 'pink', lwd = 7) # add smoothed prediction
lines(xx, pred_y, col = 'red', lwd = 1) # add raw prediction (jumps where monotonic increases in rank occur in original data)
#compare with true relationship
with(nonnorm_data, # the rank relationship  recovers the trend but not the distribution
     lines(x = x_mu, y = y_mu, col = 'cyan4', lty  = 3, lwd = 2)
     )

# Make a function ---------------------------------------------------------
PlotSpear = function(data, #data frame
                     xname = 'x', # name of x variable as a string
                     yname = 'y', # name of y variable as a string
                     ...) # plot arguments used when plotting x and y
{
  xx = with(data, 
            {
              seq(from = min(get(xname), na.rm = TRUE), 
                 to  = max(get(xname), na.rm = TRUE), 
                 length.out = 1e3)
              } )
  #Spearman rank correlation finds the correlation
  print(with(data, cor.test(get(xname), get(yname), method = 'spearman')))
  #N.B. Spearman correlation is equivalent to
  #Construct empirical cumulative density functions for x and y
  CalcECDF = function(x){efun = ecdf(x); return(efun(x))}
  data = within(data, 
                        {
                          x_cdf = CalcECDF(get(xname)) # find predicted position in range of x
                          y_cdf = CalcECDF(get(yname)) # find predicted position in range of x
                        }
  )
  #fit a linear model to predicted position (cumulative position in rank between 0 and 1)
  lm_rxy = with(data, lm(y_cdf ~ x_cdf))
  
  #generate x sequence for predictions
  #predict cumulative y value
  Ecdf_x = with(data, ecdf(get(xname)))
  pred_rxy = with(data,
                  predict(lm_rxy,
                          newdata = data.frame(x_cdf = Ecdf_x(xx)) 
                          ) 
                  )
  #predicted y value based on minimum and range of y
  pred_y = with(data, 
                {
                  diff( range(get(yname), na.rm = TRUE) )* # scale by range
                    (pred_rxy-0.5) + # subtract expected median 
                    median( get(yname), na.rm = TRUE ) # add observed median
                }
  ) 
  #smoothed prediction
  spl_y = with(data,
               {
               predict(
                  smooth.spline(x = xx, 
                                y = pred_y, 
                                nknots = round(sqrt(
                                                    sum(!is.na(get(xname)))
                                                    ))
                                ) # use a small number of knots
                        )$y # keep only y
               }
              )
  
  #plot prediction based on ranking (CDF)
  with(data, 
       {
         plot(x = get(xname), 
              y = get(yname),
              xlab = xname,
              ylab = yname,
              ...
              )
        abline(h = c(range(get(yname), na.rm = TRUE),
                     median(get(yname), na.rm = TRUE)
                     ), 
               lty = c(1,1,2)
               ) # include a dashed line of no correlation
       }
  )
  lines(xx, spl_y, col = 'pink', lwd = 7) # smoothed prediction
  lines(xx, pred_y, col = 'red', lwd = 1) # original prediction (jumps where monotonic increases in rank occur in original data)
}


# Select and load dataset of interest -------------------------------------

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


#dataset updated to include estimates of age range where not known
path_file = FT_select_file(file_type = 'dataR2.csv')

df = read.table(file = path_file,
                sep = ',', # ',' for '0.0' decimals csv, ';' for '0,0' 
                header = TRUE
)

# . Generate plots for each treatment -------------------------------------

PlotSpear(data = subset(df, Treatment == 1),
          x = 'Dspeed',
          y = 'LI',
          main = 'Treatment 1',
          pch = 21, lwd = 2,
          col = adjustcolor('darkblue', alpha = 0.7), 
          bg = adjustcolor('white', alpha.f = 0.8) 
          )
## S = 210267, p-value = 0.001369
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##   rho 
## 0.2878109 

PlotSpear(data = subset(df, Treatment == 2),
          x = 'Dspeed',
          y = 'LI',
          main = 'Treatment 2',
          pch = 21, lwd = 2,
          col = adjustcolor('orange3', alpha = 0.7), 
          bg = adjustcolor('white', alpha.f = 0.8) 
          )
## S = 233425, p-value = 0.4011
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##   rho 
## 0.07904434 


PlotSpear(data = subset(df, Treatment == 1),
          x = 'Age',
          y = 'LI',
          main = 'Treatment 1',
          pch = 21, lwd = 2,
          col = adjustcolor('darkblue', alpha = 0.7), 
          bg = adjustcolor('white', alpha.f = 0.8) 
          )
## S = 24700, p-value = 0.006161
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##   rho 
## 0.3469165 

PlotSpear(data = subset(df, Treatment == 2),
          x = 'Age',
          y = 'LI',
          main = 'Treatment 2',
          pch = 21, lwd = 2,
          col = adjustcolor('orange3', alpha = 0.7), 
          bg = adjustcolor('white', alpha.f = 0.8) 
          )

## S = 42470, p-value = 0.3644
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##   rho 
## 0.1134439 

