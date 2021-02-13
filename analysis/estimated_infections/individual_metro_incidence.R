library(epidemia)
library(rstanarm)
library(outbreakinfo)
data("EuropeCovid")
options(mc.cores = 30, warn = 1)
library(tidyverse)


state_lockdowns <- read_csv("./state_lockdown_dates.csv") %>%
    dplyr::filter(!states %in% c("District of Columbia", "Hawaii", "Puerto Rico")) %>%
    dplyr::filter(!is.na(parsed_date))

## Get metro areas with population >=1 million
metro_cases <- getMetroByCountry()
metro_cases <- metro_cases %>%
    dplyr::filter(population >= 1000000) %>%
    as_tibble()

## Add lockdown column
metro_cases <- metro_cases %>%
    dplyr::filter(date <= as.Date("2020-05-15")) %>%
    drop_na(sub_parts) %>%
    mutate(
        date = parse_date(date),
        state_name = sub_parts %>%
            map_chr(~{
                .x %>%
                    select(state_name) %>%
                    first() %>%
                    first()
            }),
        dead_numIncrease = ifelse(dead_numIncrease < 0, 0, dead_numIncrease)
    ) %>%
    dplyr::filter(state_name %in% state_lockdowns$states) %>%
    group_by(name, state_name) %>%
    group_modify(~{
        .x %>%
            mutate(
                lockdown = ifelse(date < state_lockdowns %>% dplyr::filter(states == .y$state_name) %>% select(parsed_date) %>% first(), 0, 1),
                )
    }) %>%
    ungroup() %>%
    rename(
        deaths = dead_numIncrease,
        country = name
    )

## Get metro areas within Louisiana to get "Other Louisiana" counties
shreveport <- getLocationData("Shreveport-Bossier City, LA") %>%
    as_tibble() %>%
    mutate(
        date = parse_date(date)
    )

nola <- metro_cases %>%
    dplyr::filter(country == "New Orleans-Metairie, LA")

## Get counties in the two metros in Louisiana with sampling to get "Other Louisiana"
nola_counties <- nola %>% head(1) %>% select(sub_parts) %>% first() %>% first() %>%
    select(location_id) %>%
    mutate(fips = str_replace(location_id, "USA_US-LA_", "")) %>%
    select(fips) %>%
    as_vector()

shreveport_counties <- shreveport %>% head(1) %>% select(sub_parts) %>% first() %>% first() %>%
    select(location_id) %>%
    mutate(fips = str_replace(location_id, "USA_US-LA_", "")) %>%
    select(fips) %>%
    as_vector()

other_la <- getAdmn2ByState("Louisiana") %>%
    as_tibble() %>%
    dplyr::filter(!(iso3 %in% nola_counties) & !(iso3 %in% shreveport_counties))

other_la_pop <- other_la %>%
    distinct(name, population) %>%
    summarise(n = sum(population)) %>%
    first()

## Aggregate over all other counties
other_la <- other_la %>%
    group_by(
        date
    ) %>%
    summarise(
        population = other_la_pop,
        dead_numIncrease = sum(dead_numIncrease)
    ) %>%
    ungroup() %>%
    mutate(
        name = "Other_Louisiana"
    )

la_lockdown <- as.Date("2020-03-23")

## Sum up population for Other Louisiana and Shreveport
la_metro <- other_la %>%
    mutate(
        date = parse_date(date)
    ) %>%
    bind_rows(shreveport)

## Get population for all metros
pop_tmp <- la_metro %>%
    select(name, population) %>%
    drop_na(population) %>%
    distinct() %>%
    rename(
        country = name
    ) %>%
    bind_rows(
        metro_cases %>%
        distinct(country, population)
    )

## Add lockdown column
df <- la_metro %>%
    select(date, name, dead_numIncrease) %>%
    dplyr::filter(date <= as.Date("2020-05-15")) %>%
    group_by(name) %>%
    group_modify(~{
        .x %>%
            mutate(
                lockdown = ifelse(date < la_lockdown, 0, 1),
                dead_numIncrease = ifelse(dead_numIncrease < 0, 0, dead_numIncrease)
            )
    }) %>%
    ungroup() %>%
    rename(
        deaths = dead_numIncrease,
        country = name
    )

## Merge all metros
df <- metro_cases %>%
    select(country, date, deaths, lockdown) %>%
    bind_rows(df)

df$country <- factor(df$country)

## Add dates back to Jan 1st
min_dates <- df %>% group_by(country) %>% summarise(min(date))

start_date <- as.Date("2020-01-01")
## start_date <- as.Date("2020-02-11")
buffer_dates <- map2_dfr(min_dates$country, min_dates$`min(date)`, ~{
    dates <- seq.Date(start_date, .y - 1, by="day")
    data.frame(
        date = dates,
        lockdown = 0,
        deaths = 0,
        country = .x,
        population = pop_tmp %>% dplyr::filter(country == .x) %>% select(population) %>% first()
    )
})

tmp <- df %>%
    bind_rows(buffer_dates) %>%
    select(date, deaths, country, lockdown)

tmp %>%
    mutate(
        country = as.character(country)
    ) %>%
    group_by(country) %>%
    group_walk(~{
        args <- EuropeCovid
        n <- str_replace_all(.y, "[ /]", "_")
        print(n)
        args$data <- .x %>% mutate(country = as.factor(as.character(.y)))
        args$pops <- pop_tmp %>% dplyr::filter(country == .y %>% first())
        args$algorithm <- "sampling"
        args$sampling_args <- list(iter=4e3, seed=112358, control=list(adapt_delta = 0.9))
        args$rt <- epirt(
            formula = R(country,date) ~  1 + lockdown,
            prior = normal(location=0,scale=.5),
            prior_intercept = normal(location=0,scale=2),
            prior_covariance = decov(scale=0.1)
        )
        logfile <- file(paste0("metros/", n, "_2020-10-30.log"), open="wt")
        sink(logfile, type="output")
        sink(logfile, type="message")
        fit <- do.call(epim, args)
        sink(type="output")
        sink(type="message")
        close(logfile)
        saveRDS(fit, paste0("metros/", n,"_2020-10-30.rds"))
        pdf(paste0("metro_plots/Rt_infectious_", n,"_2020-10-30.pdf"), width = 20)
        par(mfrow=c(2,2))
        print(plot_rt(fit, levels = c(20,50,80,95)))
        print(plot_obs(fit, type = "deaths", levels = c(20,50,80,95)))
        print(plot_infections(fit, levels = c(20,50,80,95)))
        dev.off()
    })
