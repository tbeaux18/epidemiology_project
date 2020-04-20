################################################################

# Parses the first table from the chicago daily reports of COVID
# cases. Change the date in %Y-%m-%d format for building the
# urls to pass into the html scraping functions

################################################################

library(tidyverse)
library(magrittr)
library(xml2)

# date in "2020-04-##" format
build_data_url <- function(latest_date) {
  base_url <- "https://www.chicago.gov/city/en/sites/covid-19/home/latest-data/"
  suffix_dates <- seq(as.Date("2020-03-23"), as.Date(latest_date), by="days")
  dates_html <- sapply(suffix_dates, paste, "html", sep=".")
  complete_urls <- sapply(base_url, paste, dates_html, sep="")
  dates_url_list <- list()
  dates_url_list[[1]] <- suffix_dates
  dates_url_list[[2]] <- complete_urls
  return(dates_url_list)
}

# pass in a url to scrape 
grab_covid_data <- function(web_url) {
  webpage_url <- web_url
  webpage <- xml2::read_html(webpage_url)
  covid_data <- rvest::html_table(webpage)[[1]] %>%
    tibble::as_tibble(.name_repair = "unique") # repair the repeated columns
  return(as.data.frame(covid_data))
}

parse_chicago_df <- function(list_of_dfs, suffix_dates) {
  incidence_count <- c(rep(NA, length(list_of_dfs)))
  for (i in 1:length(list_of_dfs)){
    date_df <- chi_cov_data[[i]]
    incidence_count[i] <- date_df[1, 2]
  }
  date_case_count_df <- data.frame(suffix_dates, incidence_count, stringsAsFactors = FALSE) %>%
    mutate(
      suffix_dates = as.Date(suffix_dates, "%Y-%m-%d"), 
      incidence_count = parse_number(incidence_count)
    ) %>% set_colnames(c("date", "covid_count"))
  return(date_case_count_df)
}

# run the functions above to parse daily incidences
chicago_urls <- build_data_url("2020-04-18")

chi_cov_data <- lapply(chicago_urls[[2]], FUN=grab_covid_data)

latest_chicago_data <- parse_chicago_df(chi_cov_data, chicago_urls[[1]])

latest_chicago_data %>%
  ggplot(aes(x = date, y = covid_count)) +
  geom_line() + geom_point()



################################################################



################################################################

library(pomp)
covid19_data <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"
us_df <- read.csv(covid19_data, stringsAsFactors = F, header = T)

cookco_data <- us_df %>%
  filter(Province_State == 'Illinois' & Admin2 == 'Cook') %>%
  select(starts_with('X')) %>%
  gather(date, covid_count) %>%
  mutate(
    date = str_replace(date, 'X', '0'),
    date = as.Date(date, '%m.%d.%y'),
    day_num = difftime(date, as.Date("2020-01-24", "%Y-%m-%d"), units="days"),
    day_num = as.numeric(day_num),
    covid_count = as.numeric(covid_count)
  )

cookco_data %>%
  ggplot(aes(x = date, y = covid_count)) +
  geom_line() + geom_point()


# 2020-03-10 starting growth
cook_data <- cookco_data %>% filter(covid_count >= 1)
fit <- lm(log(covid_count) ~ day_num, data=cook_data)
slope <- coef(fit)[2]
slope.se <- coef(summary(fit))[2,2]
1*slope.se

library(pomp)
pomp(
  data=cook_data,
  times="day_num",t0=0,
  skeleton=vectorfield(
    Csnippet("
             double incidence;
             incidence = b*S*I;
             DS = -incidence;
             DI = incidence-gamma*I;")),
  rinit=Csnippet("
                 S = S_0;
                 I = I_0;"),
  paramnames=c("b","gamma","S_0","I_0"),
  statenames=c("S","I")) -> cook_likelihood

loglik.normal <- function (params) {
  x <- trajectory(cook_likelihood, params=params)
  sum(dnorm(x=obs(cook_likelihood),mean=x["I",,],
            sd=params["sigma"],log=TRUE))
}

f3 <- function (b) {
  params <- c(S_0=5150000, I_0=5000, gamma=1, b=b, sigma=1)
  loglik.normal(params)
}

b <- seq(from=0,to=0.001,by=0.00002)
ll <- sapply(b, f3)

plot(b,ll,type='l',ylab=expression(log(L)))
b.hat <- b[which.max(ll)]
abline(v=b.hat, lty=2)


################################################################



################################################################



covid19_data <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"
us_df <- read.csv(covid19_data, stringsAsFactors = F, header = T)

cookco_data <- us_df %>%
  filter(Province_State == 'Illinois' & Admin2 == 'Cook') %>%
  select(starts_with('X')) %>%
  gather(date, covid_count) %>%
  mutate(
    date = str_replace(date, 'X', '0'),
    date = as.Date(date, '%m.%d.%y'),
    day_num = difftime(date, as.Date("2020-01-23", "%Y-%m-%d"), units="days"),
    day_num = as.numeric(day_num),
    covid_count = as.integer(covid_count)
  )

# cookco_data %>%
#   ggplot(aes(x = date, y = covid_count)) +
#   geom_line() + geom_point()


# 2020-03-10 starting growth
cook_data <- cookco_data %>% filter(covid_count >= 1)


library(deSolve)

N <- 5150000
# now specify initial values for S, I and R
init <- c(S = N - cook_data$covid_count, I = cook_data$covid_count, R = 0)
Day <- 1:(length(cook_data$covid_count))

# get a good starting condition
mod <- nls(Infected ~ a*exp(b*day), 
           start = list(a = Infected[1],
                        b = log(Infected[2]/Infected[1])))

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, { dS <- -beta/N * S * I
  dI <- beta/N * S * I - gamma * I
  dR <- gamma * I
  list(c(dS, dI, dR))
  })
}

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S/N
    dI <- beta * I * S/N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- deSolve::ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[, 3]
  sum((cook_data$covid_count - fit)^2)
}

Opt <- optim(c(2*coefficients(mod)[2]/N, coefficients(mod)[2]), RSS.SIR, method = "L-BFGS-B", lower = lower, upper = upper,
             hessian = TRUE, control = list(parscale = c(1/N,1),factr = 1))

# now find the values of beta and gamma that give the
# smallest RSS, which represents the best fit to the data.
# Start with values of 0.5 for each, and constrain them to
# the interval 0 to 1.0
Opt <- optim(c(0.5, 0.5), RSS, method="L-BFGS-B", lower=c(0,0), upper=c(1,1))

# check for convergence
Opt$message



####
####
####

library(deSolve)
library(RColorBrewer)

#https://en.wikipedia.org/wiki/Timeline_of_the_2019%E2%80%9320_Wuhan_coronavirus_outbreak#Cases_Chronology_in_Mainland_China
Infected <- c(45, 62, 121, 198, 291, 440, 571, 830, 1287, 1975, 2744, 4515)
day <- 0:(length(Infected)-1)
N <- 1400000000 #pop of china

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

# use the ordinary exponential model to determine const = beta - gamma
const <- coef(mod)[2]




RSS.SIR <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  RSS <- sum((Infected - fit)^2)
  return(RSS)
}

lower = c(0, 0)
upper = c(1, 1)  ###adjust limit because different scale 1/N

### edit: get a good starting condition
mod <- nls(Infected ~ a*exp(b*day), 
           start = list(a = Infected[1],
                        b = log(Infected[2]/Infected[1])))
optimsstart <- c(2,1)*coef(mod)[2]

set.seed(12)
Opt <- optim(optimsstart, RSS.SIR, method = "L-BFGS-B", lower = lower, upper = upper,
             hessian = TRUE)
Opt

### estimated covariance matrix of coefficients
### note the large error, but also strong correlation (nearly 1)
## note scaling with estimate of sigma because we need to use Hessian of loglikelihood
sigest <- sqrt(Opt$value/(length(Infected)-1))
solve(1/(2*sigest^2)*Opt$hessian) 

####
####  using alternative parameters
####  for this we use the function SIR2
####

optimsstart <- c(coef(mod)[2],5)
lower = c(0, 1)
upper = c(1, 10^3)  ### adjust limit because we use R0 now which should be >1

set.seed(12)
Opt2 <- optim(optimsstart, RSS.SIR2, method = "L-BFGS-B",lower=lower, upper=upper,
              hessian = TRUE, control = list(maxit = 1000, 
                                             parscale = c(10^-3,1)))
Opt2

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

################################################################



################################################################


library(EpiModel)

SEIR <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Population size
    num <- s.num + e.num + i.num + r.num
    
    # Effective contact rate and FOI from a rearrangement of Beta * c * D
    ce <- R0 / i.dur
    lambda <- ce * i.num/num
    
    dS <- -lambda*s.num
    dE <- lambda*s.num - (1/e.dur)*e.num
    dI <- (1/e.dur)*e.num - (1 - cfr)*(1/i.dur)*i.num - cfr*(1/i.dur)*i.num
    dR <- (1 - cfr)*(1/i.dur)*i.num
    
    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, but within the containing list
    list(c(dS, dE, dI, dR, 
           se.flow = lambda * s.num,
           ei.flow = (1/e.dur) * e.num,
           ir.flow = (1 - cfr)*(1/i.dur) * i.num,
           d.flow = cfr*(1/i.dur)*i.num),
         num = num,
         i.prev = i.num / num,
         ei.prev = (e.num + i.num)/num)
  })
}



# ICM models EpiModel Tutorial
param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
init <- init.icm(s.num = 500, i.num = 1)
control <- control.icm(type = "SI", nsims = 10, nsteps = 300)
(mod <- icm(param, init, control))
summary(mod, at = 125)

head(as.data.frame(mod, out = "mean"))

plot(mod)

plot(mod, y = "i.num", sim.lines = TRUE, mean.smooth = FALSE, qnts.smooth = FALSE)


# Ebola model information
# reproductive number decayed logistically over time
# weekly incidence data
# We developed a discrete-time stochastic compartmental model with a time-varying reproductive number and a time step of 1 day
# nsteps = 52 weeks, 10 independent runs
control <- control.icm(type = "SIR", nsteps = 365, nsims = 10)

# figure out the simulated population of individuals
# look at illinois statistics

# define the 3 groups, can add more attributes with extension models
# chicago pop 2706000 from 2018
# cook pop 5150000 
init <- init.icm(s.num = 997, i.num = 3, r.num = 0)

# parameters, figure out how to estimate these
# act.rate = exposure-to-infection rate set to 10 times per day
# inf.prob = overall probability of infection across all those exposures set to 0.05
# rec.rate = 1/20
# dr.rate = death rate, figure out chicago's weekly death rate
# a.rate = (13.4/365)/1000) from 2017
param <- param.icm(inf.prob = 0.05, act.rate = 10, rec.rate = 1/20, 
                   a.rate = (10.5/365)/1000, ds.rate = (7/365)/1000, di.rate = (14/365)/1000, 
                   dr.rate = (7/365)/1000)

sim <- icm(param, init, control)

# death rate per 100,000 people
chicago_morb_mort <- read.csv('CDPH_Morbidity_and_Mortality.csv', stringsAsFactors = F, header = T)
