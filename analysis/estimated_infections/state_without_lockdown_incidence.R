library(epidemia)
library(rstanarm)
library(outbreakinfo)
data("EuropeCovid")
options(mc.cores = 20)
library(tidyverse)

state_lockdowns <- read_csv("./state_lockdown_dates.csv") %>%
    dplyr::filter(is.na(parsed_date))

## Get deaths
us_states <- getAdmn1ByCountry("United States of America") %>%
    dplyr::filter(name %in% state_lockdowns$states)

## Sum up population
pop_tmp <- us_states %>%
    select(name, population) %>%
    drop_na(population) %>%
    distinct() %>%
    rename(
        country = name,
        pop = population
    )

df <- us_states %>%
    select(date, name, dead_numIncrease) %>%
    mutate(date = parse_date(date)) %>%
    dplyr::filter(date <= as.Date("2020-05-15")) %>%
    group_by(name) %>%
    group_modify(~{
        .x %>%
            mutate(
                dead_numIncrease = ifelse(dead_numIncrease < 0, 0, dead_numIncrease)
            )
    }) %>%
    ungroup() %>%
    rename(
        deaths = dead_numIncrease,
        country = name
    )

df$country <- factor(df$country)

## Add dates back to Jan 1st
min_dates <- df %>% group_by(country) %>% summarise(min(date))

start_date <- as.Date("2020-01-01")
buffer_dates <- map2_dfr(min_dates$country, min_dates$`min(date)`, ~{
    dates <- seq.Date(start_date, .y - 1, by="day")
    data.frame(
        date = dates,
        deaths = 0,
        country = .x,
        population = pop_tmp %>% dplyr::filter(country == .x) %>% select(pop) %>% first()
    )
})

tmp <- df %>%
    bind_rows(buffer_dates) %>%
    select(date, deaths, country)

# collect arguments for 'epim'tmp
args <- EuropeCovid
print(names(args))

args$data <- tmp
args$pops <- pop_tmp

# use NUTs sampler to fit the model
args$algorithm <- "sampling"
args$sampling_args <- list(iter=10000,seed=112358, chains = 4)

# model for the reproduction number
args$rt <- epirt(
  formula = R(country,date) ~ (1 | country),
  prior = normal(location=0,scale=.5),
  prior_intercept = normal(location=0,scale=2)
)

args$group_subset <- state_lockdowns$states
## args$group_subset <- c("California", "Illinois")
fit <- do.call(epim, args)
saveRDS(fit, "us_states_without_lockdown_2020-10-07.rds")

## pdf("../Rt_infection_us_states_2020-09-19.pdf", width = 20)
## par(mfrow=c(2,2))
## plot_rt(fit)
## plot_infections(fit)
## plot_obs(fit, type = "deaths")
## dev.off()

