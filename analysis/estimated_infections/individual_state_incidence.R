library(epidemia)
library(rstanarm)
library(outbreakinfo)
data("EuropeCovid")
options(mc.cores = 30, warn = 1)
library(tidyverse)


state_lockdowns <- read_csv("./state_lockdown_dates.csv") %>%
    dplyr::filter(!states %in% c("District of Columbia", "Hawaii", "Puerto Rico")) %>%
    dplyr::filter(!is.na(parsed_date))

## Get deaths
us_states <- getAdmn1ByCountry("United States of America") %>%
    dplyr::filter(name %in% state_lockdowns$states & name == "California")
    ## %>%
    ## dplyr::filter(name %in% c("California", "New York", "Florida", "Illinois", "Texas", "Washington", "Louisiana", "Connecticut", "Idaho"))

## data.frame(
##     names = c("California", "New York", "Florida", "Illinois", "Texas", "Washington", "Louisiana", "Connecticut", "Idaho"),
##     dates = as.Date(c("2020-03-19", "2020-03-22", "2020-04-03", "2020-03-21", "2020-04-02", "2020-03-23", "2020-03-23", "2020-03-23", "2020-03-25"))
## )

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
                lockdown = ifelse(date < state_lockdowns %>% dplyr::filter(states == .y$name) %>% select(parsed_date) %>% first(), 0, 1),
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

start_date <- as.Date("2019-12-01")
buffer_dates <- map2_dfr(min_dates$country, min_dates$`min(date)`, ~{
    dates <- seq.Date(start_date, .y - 1, by="day")
    data.frame(
        date = dates,
        lockdown = 0,
        deaths = 0,
        country = .x,
        population = pop_tmp %>% dplyr::filter(country == .x) %>% select(pop) %>% first()
    )
})

tmp <- df %>%
    bind_rows(buffer_dates) %>%
    select(date, deaths, country, lockdown)

.x <- tmp %>%
    mutate(
        country = as.character(country)
    ) %>%
    dplyr::filter(country == "California")

tmp %>%
    mutate(
        country = as.character(country)
    ) %>%
    dplyr::filter(country == "Florida") %>%
    group_by(country) %>%
    group_walk(~{
        args <- EuropeCovid
        n <- .y %>% first() %>% str_replace_all("[ /,]", "_")
        print(n)
        args$data <- .x %>% mutate(country = as.factor(as.character(.y)))
        args$pops <- pop_tmp %>% dplyr::filter(country == .y %>% first())
        args$algorithm <- "sampling"
        args$sampling_args <- list(iter=2e3, seed=112358, control=list(adapt_delta = 0.9))
        args$rt <- epirt(
            formula = R(country,date) ~  1 + lockdown,
            prior = normal(location=0,scale=.5),
            prior_intercept = normal(location=0,scale=2),
            prior_covariance = decov(scale=0.1)
        )
        logfile <- file(paste0("states/",n, "_2020-10-30.log"), open="wt")
        sink(logfile, type="output")
        sink(logfile, type="message")
        fit <- do.call(epim, args)
        sink(type="output")
        sink(type="message")
        close(logfile)
        saveRDS(fit, paste0("states/",n,"_2020-10-30.rds"))
        pdf(paste0("state_plots/Rt_infectious_", n,"_2020-10-30.pdf"), width = 20)
        par(mfrow=c(2,2))
        print(plot_rt(fit, levels = c(20,50,80,95)))
        print(plot_obs(fit, type = "deaths", levels = c(20,50,80,95)))
        print(plot_infections(fit, levels = c(20,50,80,95)))
        dev.off()
    })
