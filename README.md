Multivariate COM-Poisson distribution: Sarmanov method and Bayesian
inference
================

Before installing `multcp` ensure that C++ compiler is properly set up
to be used with R. We suggest to test this with
`install.packages("RcppArmadillo") library(RcppArmadillo)` as our
package links to this framework. Once this is done, `multcp` can be
installed from Git Hub with the aid of `devtools` and is ready to be
used.

``` r
library(devtools)
install_github("luizapiancastelli/multcp")
library(multcp)
```

### Simulating MCOMP data

Data can be simulated from the MCOMP sampler algorithm in Section 4.4
with the `rmcomp` function that takes in the arguments: the number of
random draws `n`, parameter vectors `lambda` and `nu` of lenght `d`,
`d`x`d` matrix `delta` and scalar `omega`. `N_r` is the number of
auxiliary draws used to estimate the intractable ratios via importance
sampling (Section 4.1). It is also possible to supply the truncation
tolerance for computing the coniditional probability functions ε, here
named `tol`.

``` r
n = 100
lambda = c(1.5, 1, 0.5)
nu = c(1, 0.5, 2)
omega = 1
delta = matrix(0, ncol = length(lambda), nrow = length(lambda)) #Other than upper triangle elements should be 0
delta[1,2] = 2.5
delta[1,3] = 2
delta[2,3] = 3
N_r = 10000

Y = rmcomp(n, lambda, nu, delta, omega, N_r)
colMeans(Y)
```

    ## [1] 1.51 1.61 0.27

``` r
cor(Y)
```

    ##           [,1]      [,2]      [,3]
    ## [1,] 1.0000000 0.2086069 0.1177070
    ## [2,] 0.2086069 1.0000000 0.3019776
    ## [3,] 0.1177070 0.3019776 1.0000000

### GIMH

The pseudo-marginal approach (Section 5.1) is implemented in function
`GIMH` that runs (possibly in parallel with `ncores`) `chains` of
`n_iter` iterations and `burn_in` period. Estimation of intractable
constants is done with `N_r` (integer) draws for the ratios, and `N_z`
can be either supplied as integer or calibrated. Calibration is done
targetting the desired log-likelihood standard deviation `target_sd`,
starting from increasing the initial number `init` of draws until this
is met. It is recommended to set `target_sd` between 1.2-1.7 for good
mixing, but more details con be found in Section 5.1.1. Prior parameters
are supplied via a `priors_list`, and parallelization can be disabled by
setting `ncores=1`. The following example uses the simulated data set of
configuration 1 (bivariate with positive correlation), with the other
settings available named `Y_c2, Y_c3, Y_c4`.

``` r
data("Y_c1") #Loads the simulated data of configuration 1

N_aux_r = 10000
N_aux_z = list('target_sd' =1.2, 'init' = 10000) #Will run calibration of Nz
burn_in = 10000
n_iter = 2000

priors_list = list('lambda' = list('sd' = 100),
                   'nu' = list('sd' = 100),
                   'omega' = list('sd' = 5))

#Progress is printed to the console if single chain, for multiple chains refer to a .txt file that is saved
#to the working directory every 100 iterations
run_gimh = GIMH(Y_c1, burn_in, n_iter, N_aux_r, N_aux_z, priors_list, chains = 2, ncores = 2)

names(run_gimh) #Returns a list with timing and raw MCMC object
```

The function returns a list with the total algorithm’s runtime `$time`
and a raw MCMC object `$mcmc` of length `chains`. The former contains
the posterior draws, acceptance rates, initial values and proposal
parameters and can be processed and summarised using `process_mcmc` and
`posterior_summaries`.

``` r
names(run_gimh$mcmc[[1]]) #Raw output for each of the {1,...,chains}

process_gimh = process_mcmc(run_gimh$mcmc)
gimh_summary = posterior_summaries(run_gimh$mcmc)

gimh_summary$Rhat     #Rhat statistic
gimh_summary$Stats    #Some posterior summaries using the combined MCMC draws

gimh_summary$Density_plots  #Posterior density plots
gimh_summary$Trace_plots    #Trace plots
```

### Noisy Exchange algorithm

The `Exchange` function works similarly, providing the methodology
introduced in Section 5.2. Here the value of ε can be supplied as `N_z`
no longer applies. The output is in the same format as `GIMH` and the
functions `process_mcmc` and `posterior_summaries` can also be
used.

``` r
run_exchange = Exchange(n_iter, burn_in, Y_c1, N_aux_r, priors_list, chains=2,  ncores = 2, tol = 0.001)

exchange_summary = posterior_summaries(run_exchange$mcmc)
exchange_summary$Stats
```

### Premier League Regression

Our implementation of the Exchange algorithm for the regression model of
Section 4 is available from the function `Exchange_regression` and the
Premier League data is loaded with `data("Y_premier")`.

``` r
data("Y_premier")
Y = as.matrix(Y_premier[,1:2])
X =Y_premier[,3]

burn_in = 10000
n_iter = 2000
N_aux_r = 10000

priors_list = list('gamma0' = list('sd' = 100), 
                   'gamma' = list('sd' = 100),
                   'nu' = list('sd' = 100),
                   'omega' = list('sd' = 5))

run_premier_exchange = Exchange_regression(burn_in, n_iter, Y, X, N_aux_r, priors_list, chains=2,  ncores = 2, tol = 0.001)

application_summary = posterior_summaries(run_premier_exchange$mcmc) #Summaries of MCMC output
application_summary$Stats
application_summary$Density_plots
application_summary$Trace_plots
```
