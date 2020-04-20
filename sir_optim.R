library(deSolve)
library(RColorBrewer)


covid19_data <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"
us_df <- read.csv(covid19_data, stringsAsFactors = F, header = T)

cook_data <- us_df %>%
  filter(Province_State == 'Illinois' & Admin2 == 'Cook') %>%
  dplyr::select(starts_with('X')) %>%
  gather(date, covid_count) %>%
  mutate(
    date = str_replace(date, 'X', '0'),
    date = as.Date(date, '%m.%d.%y'),
    day_num = difftime(date, as.Date("2020-01-23", "%Y-%m-%d"), units="days"),
    day_num = as.numeric(day_num),
    covid_count = as.integer(covid_count)
  ) %>% filter(covid_count >= 4)

N <- 5150000
# now specify initial values for S, I and R
Infected <- cook_data %>% pull(covid_count)

day <- 0:(length(Infected)-1)
init <- c(S = N-Infected[1], I = Infected[1], R = 0)

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

### plotting different values R0

const <- modnested$const
R0 <- modnested$R0

# graph
plot(-100,-100, xlim=c(0,80), ylim = c(1,6*10^4), log="", 
     ylab = "infected", xlab = "days")
title(bquote(paste("scenario's for different ", R[0])), cex.main = 1)

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

