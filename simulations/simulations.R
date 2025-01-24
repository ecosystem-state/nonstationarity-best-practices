library(MARSS)
library(marssTMB)

set.seed(123)
pars <- expand.grid(seed = 1:200,
                    n_years = c(25,50,100, 200),
                    covar_trend = c(0, 0.01, 0.02),
                    coef_sd = c(0.01,0.1,1), # variability of random walk of the time varying covariate
                    rho = c(-0.9,-0.5,0,0.5,0.9)) # autocorrelation of cov time series. assume to have sd =1
pars$R <- NA
pars$q_alpha <- NA
pars$q_beta <- NA
pars$convergence <- NA
pars$rho_cov <- NA

for(i in 1:nrow(pars)) {
  
set.seed(pars$seed[i])
## number of years of data
TT <- pars$n_years[i]
years <- 1:TT
coef_sd <- pars$coef_sd[i]

## get predictor variable and z-score

tv_coef <- matrix(cumsum(rnorm(TT, mean = 0, sd = coef_sd)), nrow=1) # time varying effect
#tv_coef_z <- matrix((tv_coef - mean(tv_coef))/sd(tv_coef), nrow = 1)

x0 <- cumsum(rnorm(TT, mean = 0, sd = 0.1))

x <- c(arima.sim(n = TT, list(ar = pars$rho[i]), sd = 1))

## get response variable: logit(survival)
dat <- x0 + tv_coef * x + seq(0,TT-1)*pars$covar_trend[i] + rnorm(TT, 0, 0.01) 

## number of regr params (slope + intercept)
m <- dim(x_z)[1] + 1

## for process eqn
B <- diag(m)  ## 2x2; Identity
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now
diag(Q) <- c("q.alpha", "q.beta")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- x  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0
R <- matrix("r")  ## 1x1; scalar = r

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

dlm_1 <- MARSS(dat, inits = inits_list, model = mod_list, method="TMB")

#tidy_pars <- tidy(dlm_1)
pars$R[i] <- dlm_1$par$R
pars$q_alpha[i] <- sqrt(dlm_1$par$Q[1,1])
pars$q_beta[i] <- sqrt(dlm_1$par$Q[2,1])
pars$convergence[i] <- dlm_1$convergence

pars$rho_cov[i] <- cor(dlm_1$states[2,], c(tv_coef))
print(i)
}

saveRDS(pars, "simulation_pars.rds")
