

library(tidyverse)
library(magrittr)
library(incidence)


us_data_il <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"

il_df <- read.csv(us_data_il, stringsAsFactors = F, header = T)

illinois_data <- il_df %>%
  filter(Province_State == 'Illinois' & Admin2 == 'Cook') %>%
  dplyr::select(starts_with('X')) %>%
  gather(date, cases) %>%
  filter(cases > 2) %>%
  mutate(
    date = str_replace(date, 'X', '0'),
    date = as.Date(date, '%m.%d.%y')
  )
illinois_data %<>% mutate(day_num = seq(1, nrow(illinois_data), 1.0))

dates <- rep(illinois_data$date, times=illinois_data$cases)
i <- incidence(dates)
f <- fit(i)

model <- f$model
plot(i, show_cases = TRUE, fit = f)

newdata <- data.frame(dates.x=seq(51.5, 81.5, 1.0))

x_dates <- seq(as.Date("2020-04-21"), length.out = 31, by="day")

pred_df <- data.frame(predict(model, newdata, type='response', interval='confidence'), dates.x=newdata, type='prediction', dates=x_dates)
model_df <- data.frame(fitted=model$fitted.values, dates.x=model$model$dates.x, obs.vals=model$model$`log(counts)`, date=illinois_data$date)


library(projections)
library(distcrete)
si <- distcrete("gamma", interval = 1L, shape = 3.96, w = 0.5)
p <- project(i, R = 1.08, si = si, n_sim = 10, n_days = 30, R_fix_within = TRUE)

pred_count <- as.data.frame(p) %>%
  gather(sim, count, -dates) %>%
  group_by(dates) %>%
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n()
  ) %>%
  mutate(
    type = 'predicted',
    se_count = sd_count / sqrt(n),
    lower_ci = mean_count - qt(1 - (0.05 / 2), n - 1) * se_count,
    upper_ci = mean_count + qt(1 - (0.05 / 2), n - 1) * se_count
  )

pred_count %>%
  ggplot(aes(x=dates, y=mean_count)) +
  geom_line(color='blue', size=1) +
  geom_line(data=illinois_data, aes(x=date, y=cases), color='black', size=1) +
  geom_line(data=model_df, aes(x=date,y=exp(fitted)), color='red', size=1) +
  scale_x_date(date_labels = "%b/%d", date_breaks="2 week") +
  scale_y_continuous(breaks=seq(0,150000,10000)) +
  labs(x = 'Time', y='Daily Incidence', title='Predicted Incidence Count') +
  theme(plot.title = element_text(hjust = 0.5))




########################################################################

# Estimating reproductive number from incidence data


########################################################################

# r_t_range is a vector of possible values for R_t
R_T_MAX = 12
r_t_range = seq(0, R_T_MAX, length = R_T_MAX*100 + 1)

# Gamma is 1/serial interval
# https://wwwnc.cdc.gov/eid/article/26/6/20-0357_article
GAMMA = 1/4

library(HDInterval)
library(smoother)
#' Compute new cases and smooth them
smooth_new_cases <- function(cases){
  cases %>%
    arrange(date) %>%
    mutate(new_cases = c(cases[1], diff(cases))) %>%
    mutate(new_cases_smooth = round(
      smoother::smth(new_cases, window = 7, tails = TRUE)
    )) %>%
    select(date, new_cases, new_cases_smooth)
}

compute_likelihood <- function(cases){
  likelihood <- cases %>%
    filter(new_cases_smooth > 0) %>%
    mutate(
      r_t = list(r_t_range),
      lambda = map(lag(new_cases_smooth, 1), ~ .x * exp(GAMMA * (r_t_range - 1))),
      likelihood_r_t = map2(new_cases_smooth, lambda, dpois, log = TRUE)
    ) %>%
    slice(-1) %>%
    select(-lambda) %>%
    unnest(c(likelihood_r_t, r_t))
}

compute_posterior <- function(likelihood){
  likelihood %>%
    arrange(date) %>%
    group_by(r_t) %>%
    mutate(posterior = exp(
      zoo::rollapplyr(likelihood_r_t, 7, sum, partial = TRUE)
    )) %>%
    group_by(date) %>%
    mutate(posterior = posterior / sum(posterior, na.rm = TRUE)) %>%
    # HACK: NaNs in the posterior create issues later on. So we remove them.
    mutate(posterior = ifelse(is.nan(posterior), 0, posterior)) %>%
    ungroup() %>%
    select(-likelihood_r_t)
}

estimate_rt <- function(posteriors){
  posteriors %>%
    group_by(date) %>%
    summarize(
      r_t_simulated = list(sample(r_t_range, 10000, replace = TRUE, prob = posterior)),
      r_t_most_likely = r_t_range[which.max(posterior)]
    ) %>%
    mutate(
      r_t_lo = map_dbl(r_t_simulated, ~ hdi(.x)[1]),
      r_t_hi = map_dbl(r_t_simulated, ~ hdi(.x)[2])
    ) %>%
    select(-r_t_simulated)
}

plot_estimates <- function(estimates){
  estimates %>%
    ggplot(aes(x = date, y = r_t_most_likely)) +
    geom_point(color='darkorange', alpha=0.9, size=3) +
    geom_line(color = 'steelblue', size=2, alpha=0.8) +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    geom_ribbon(aes(ymin = r_t_lo, ymax = r_t_hi), fill = 'darkred', alpha = 0.2) +
    labs(title = 'R(t) Estimates', x = 'Time', y = 'R(t) estimate', subtitle = 'Illinois') +
    scale_x_date(date_labels = "%b/%d", date_breaks="1 week") +
    coord_cartesian(ylim = c(0, 4)) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
}

il_estimates <- illinois_data %>%
  smooth_new_cases() %>%
  compute_likelihood() %>%
  compute_posterior() %>%
  estimate_rt()

il_estimates %>% plot_estimates()


########################################################################

# Estimating beta and gamma parameters from current data

########################################################################

library(deSolve)
library(RColorBrewer)



# now specify initial values for S, I and R
Infected <- illinois_data %>% pull(cases)
day <- 0:(length(Infected)-1)
N <- 5000000 #pop

###edit 1: use different boundary condiotion
###init <- c(S = N-1, I = 1, R = 0)
init <- c(S = N-Infected[1], I = Infected[1], R = 0)
plot(day, Infected)

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  ####edit 2; use equally scaled variables 
  with(par, { dS <- -beta * (S/N) * I
  dI <- beta * (S/N) * I - gamma * I
  dR <- gamma * I
  list(c(dS, dI, dR))
  })
}

SIR2 <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  ####
  #### use as change of variables variable
  #### const = (beta-gamma)
  #### delta = gamma/beta
  #### R0 = beta/gamma > 1 
  #### 
  #### beta-gamma = beta*(1-delta)
  #### beta-gamma = beta*(1-1/R0)
  #### gamma = beta/R0
  with(par, { 
    beta  <- const/(1-1/R0)  
    gamma <- const/(R0-1)  
    dS <- -(beta * (S/N)      ) * I 
    dI <-  (beta * (S/N)-gamma) * I 
    dR <-  (             gamma) * I
    list(c(dS, dI, dR))
  })
}

RSS.SIR2 <- function(parameters) {
  names(parameters) <- c("const", "R0")
  out <- ode(y = init, times = day, func = SIR2, parms = parameters)
  fit <- out[ , 3]
  RSS <- sum((Infected - fit)^2)
  return(RSS)
}

### plotting different values R0

library(minpack.lm)
mod <- nlsLM(Infected ~ a*exp(b*day), 
             start = list(a = Infected[1],
                          b = log(Infected[2]/Infected[1])))

optimsstart <- c(2,1)*coef(mod)[2]

# use the ordinary exponential model to determine const = beta - gamma
const <- coef(mod)[2]

RSS.SIR <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  RSS <- sum((Infected - fit)^2)
  return(RSS)
}

optimsstart <- c(coef(mod)[2],5)
lower = c(0, 1)
upper = c(1, 10^3)  ### adjust limit because we use R0 now which should be >1

set.seed(12)
Opt2 <- optim(optimsstart, RSS.SIR2, method = "L-BFGS-B",lower=lower, upper=upper,
              hessian = TRUE, control = list(maxit = 1000, 
                                             parscale = c(10^-3,1)))
Opt2_par <- setNames(Opt2$par, c("beta", "gamma"))

t <- 1:70
fit <- data.frame(ode(y = init, times = t, func = SIR2, parms = Opt2_par))
col <- 1:3
matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 2, lty = 1, col = col, log = "y")
points(day, log(Infected))
legend("bottomright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 2, col = col, inset = 0.05)
# now the estimated variance of the 1st parameter is small
# the 2nd parameter is still with large variance
#
# thus we can predict beta - gamma very well
# this beta - gamma is the initial growth coefficient
# but the individual values of beta and gamma are not very well known
#
# also note that hessian is not at the MLE since we hit the lower boundary
#
sigest <- sqrt(Opt2$value/(length(Infected)-1))
solve(1/(2*sigest^2)*Opt2$hessian)

#### We can also estimated variance by
#### Monte Carlo estimation
##
## assuming data to be distributed as mean +/- q mean
## with q such that mean RSS = 52030
##
## 
##


### Two functions RSS to do the optimization in a nested way
RSS.SIRMC2 <- function(const,R0) {
  parameters <- c(const=const, R0=R0)
  out <- ode(y = init, times = day, func = SIR2, parms = parameters)
  fit <- out[ , 3]
  RSS <- sum((Infected_MC - fit)^2)
  return(RSS)
}
RSS.SIRMC <- function(const) {
  optimize(RSS.SIRMC2, lower=1,upper=10^5,const=const)$objective
}

getOptim <- function() {
  opt1 <- optimize(RSS.SIRMC,lower=0,upper=1)
  opt2 <- optimize(RSS.SIRMC2, lower=1,upper=10^5,const=opt1$minimum)
  return(list(RSS=opt2$objective,const=opt1$minimum,R0=opt2$minimum))
}

# modeled data that we use to repeatedly generate data with noise
Opt_par <- Opt2$par
names(Opt_par) <- c("const", "R0")
modInfected <- data.frame(ode(y = init, times = day, func = SIR2, parms = Opt_par))$I


modInfected_df <- data.frame(ode(y = init, times = day, func = SIR2, parms = Opt_par))



# doing the nested model to get RSS
set.seed(1)
Infected_MC <- Infected
modnested <- getOptim()

errrate <- modnested$RSS/sum(Infected) 


par <- c(0,0)
for (i in 1:100) {
  Infected_MC <- rnorm(length(modInfected),modInfected,(modInfected*errrate)^0.5)
  OptMC <- getOptim()
  par <- rbind(par,c(OptMC$const,OptMC$R0))
}
par <- par[-1,]

plot(par, xlab = "const",ylab="R0",ylim=c(1,1))
title("Monte Carlo simulation")
cov(par)


###conclusion: the parameter R0 can not be reliably estimated

##### End of Monte Carlo estimation


### plotting different values R0

# use the ordinary exponential model to determine const = beta - gamma
const <- coef(mod)[2]
R0 <- 1.1

# graph
plot(-100,-100, xlim=c(0,80), ylim = c(1,N), log="y", 
     ylab = "infected", xlab = "days", yaxt = "n")
axis(2, las=2, at=10^c(0:9),
     labels=c(expression(1),
              expression(10^1),
              expression(10^2),
              expression(10^3),
              expression(10^4),
              expression(10^5),
              expression(10^6),
              expression(10^7),
              expression(10^8),
              expression(10^9)))
axis(2, at=rep(c(2:9),9)*rep(10^c(0:8),each=8), labels=rep("",8*9),tck=-0.02)
title(bquote(paste("scenario's for different ", R[0])), cex.main = 1)

# time
t <- seq(0,60,0.1)

# plot model with different R0
for (R0 in c(1.1,1.2,1.5,2,3,5,10)) {
  fit <- data.frame(ode(y = init, times = t, func = SIR2, parms = c(const,R0)))$I
  lines(t,fit)
  text(t[601],fit[601],
       bquote(paste(R[0], " = ",.(R0))),
       cex=0.7,pos=4)
}

# plot observations
points(day,Infected)


# model function
SIR2 <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, { 
    beta  <- const/(1-1/R0)  
    gamma <- const/(R0-1)  
    dS <- -(beta * (S/N)      ) * I 
    dI <-  (beta * (S/N)-gamma) * I 
    dR <-  (             gamma) * I
    list(c(dS, dI, dR))
  })
}

### Two functions RSS to do the optimization in a nested way
RSS.SIRMC2 <- function(R0,const) {
  parameters <- c(const=const, R0=R0)
  out <- ode(y = init, times = day, func = SIR2, parms = parameters)
  fit <- out[ , 3]
  RSS <- sum((Infected_MC - fit)^2)
  return(RSS)
}

RSS.SIRMC <- function(const) {
  optimize(RSS.SIRMC2, lower=1,upper=10^5,const=const)$objective
}

# wrapper to optimize and return estimated values
getOptim <- function() {
  opt1 <- optimize(RSS.SIRMC,lower=0,upper=1)
  opt2 <- optimize(RSS.SIRMC2, lower=1,upper=10^5,const=opt1$minimum)
  return(list(RSS=opt2$objective,const=opt1$minimum,R0=opt2$minimum))
}

# doing the nested model to get RSS
Infected_MC <- Infected
modnested <- getOptim()

rss <- sapply(seq(0.3,0.5,0.01), 
              FUN = function(x) optimize(RSS.SIRMC2, lower=1,upper=10^5,const=x)$objective)

plot(seq(0.3,0.5,0.01),rss)

optimize(RSS.SIRMC2, lower=1,upper=10^5,const=0.35)

# view
modnested


# plot model with different R0
t <- seq(0,50,0.1)
for (R0 in c(modnested$R0,1.07,1.08,1.09,1.1,1.11)) {
  fit <- data.frame(ode(y = init, times = t, func = SIR2, parms = c(const,R0)))$I
  lines(t,fit,col=1+(modnested$R0==R0))
  text(t[501],fit[501],
       bquote(paste(R[0], " = ",.(R0))),
       cex=0.7,pos=4,col=1+(modnested$R0==R0))
}

# plot observations
points(day,Infected, cex = 0.7)

########################################################################

# Simulating COVID19 outbreak based on incidence data


########################################################################

library(EpiModel)

l_rt <-  il_estimates$r_t_most_likely[nrow(il_estimates)]
# param <- param.dcm(R0 = l_rt, e.dur = 10, i.dur = 14, cfr = c(0.5, 0.7, 0.9))
# 
# init <- init.dcm(s.num = 5150000, e.num = 1000000, i.num = 0, r.num = 0,
#                  se.flow = 0, ei.flow = 0, ir.flow = 0, d.flow = 0)
# control <- control.dcm(nsteps = 100, dt = 1, new.mod = SEIR)
# mod <- dcm(param, init, control)
# plot(mod, y = "si.flow", lwd = 4, col = "firebrick",
#      main = "Disease Incidence", legend = "n")

# inf.prob and rec.rate 
# act rate 
param <- param.dcm(inf.prob = 0.2, act.rate = 1.3, rec.rate = 0.2)
init <- init.dcm(s.num = 5150000, i.num = 5000, r.num = 0)
control <- control.dcm(type = "SIR", nsteps = 100, dt = 1, nsims=10)
mod <- dcm(param, init, control)

sim_df <- as.data.frame(mod)

sim_df %<>% mutate(date = seq(as.Date("2020-03-01"), as.Date("2020-06-08"), 1))

sim_df %>%
  ggplot(aes(x = date, y = si.flow)) +
  geom_line(color='blue', size=1.5) +
  geom_line(data=sim_df, aes(x = date, y=i.num), color='black', size=1.5)+
  geom_line(data=illinois_data, aes(x=date, y=cases), color='red', size=1.5) +
  scale_x_date(date_labels = "%b/%d", date_breaks="2 week") +
  scale_y_continuous(breaks=seq(0,200000,10000), labels = function(x) format(x, scientific = TRUE)) +
  labs(x = 'Time', y = 'Number Infected', title='Simulated Prediction Model')

