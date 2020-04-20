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

library(EpiModel)

fit1 <- lm(log(measles)~biweek, data=latest_chicago_data)
summary(fit1)





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
control <- control.icm(type = "SIR", nsteps = 52, nsims = 10)

# figure out the simulated population of individuals
# look at illinois statistics

# define the 3 groups, can add more attributes with extension models
init <- init.icm(s.num = 997, i.num = 3, r.num = 0)

# parameters, figure out how to estimate these
# act.rate = exposure-to-infection rate set to 10 times per day
# inf.prob = overall probability of infection across all those exposures set to 0.05
# rec.rate = 1/20
# dr.rate = death rate, figure out chicago's weekly death rate

param <- param.icm(inf.prob = 0.05, act.rate = 10, rec.rate = 1/20, 
                   a.rate = (10.5/365)/1000, ds.rate = (7/365)/1000, di.rate = (14/365)/1000, 
                   dr.rate = (7/365)/1000)

sim <- icm(param, init, control)