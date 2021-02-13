library(tidyverse)
library(lubridate)
library(outbreakinfo)
library(epidemia)

state_lockdowns <- read_csv("./state_lockdown_dates.csv") %>%
    dplyr::filter(!states %in% c("District of Columbia", "Hawaii", "Puerto Rico"))

## Get reported cases
cases <- getAdmn1ByCountry("United States of America") %>%
    dplyr::filter(name %in% state_lockdowns$states) %>%
    mutate(
        date = parse_datetime(date)
    ) %>%
    as_tibble()

## Get reported for two metro areas
la_metros <- c("New Orleans-Metairie, LA", "Shreveport-Bossier City, LA")
metro_cases <- getLocationData(la_metros) %>%
    mutate(
        date = parse_datetime(date)
    ) %>%
    as_tibble()

## Get counties in the two metros to get Other Louisiana
nola_shreveport_counties <- metro_cases %>%
    distinct(sub_parts) %>%
    map_dfr(~{.x}) %>%
    select(location_id) %>%
    mutate(
        location_id = str_replace(location_id, "USA_US-LA_", "")
    ) %>%
    as_vector()

other_la_cases <- getAdmn2ByState("Louisiana") %>%
    as_tibble() %>%
    dplyr::filter(!(iso3 %in% nola_shreveport_counties))

other_la_cases_pop <- other_la_cases %>%
    distinct(name, population) %>%
    summarise(n = sum(population)) %>%
    first()

## Aggregate over all other counties
other_la_cases <- other_la_cases %>%
    group_by(
        date
    ) %>%
    summarise(
        population = other_la_cases_pop,
        confirmed = sum(confirmed, na.rm=TRUE),
        confirmed_numIncrease = sum(confirmed_numIncrease, na.rm = TRUE),
        dead_numIncrease = sum(dead_numIncrease, na.rm=TRUE),
        dead = sum(dead, na.rm = TRUE),
    ) %>%
    ungroup() %>%
    mutate(
        name = "Other_Louisiana",
        date = parse_datetime(date)
    )

## Merge all
cases <- bind_rows(cases, metro_cases, other_la_cases)

## Get all modelled infections for states with and without lockdown and LA metro areas
rds_paths <- list.files(c("./states", "./state_without_lockdown"), pattern="*.rds", full.names=TRUE)
rds_paths <- c("./metros/Shreveport-Bossier_City,_LA_2020-10-30.rds", "./metros/New_Orleans-Metairie,_LA_2020-10-30.rds", "./metros/Other_Louisiana_2020-10-30.rds", rds_paths) #Add LA metro areas
infections <- rds_paths  %>% map_dfr(~{
    fit <- readRDS(.x)
    inf <- posterior_infections(fit)
    tibble(
        time = inf$time,
        group = inf$group,
        median = apply(inf$draws, 2, median),
        lower = apply(inf$draws, 2, function(x){
            quantile(x, 0.025)
        }),
        upper = apply(inf$draws, 2, function(x){
            quantile(x, 0.975)
        })
    )
})

## Plot infections
states_without_lockdown <- state_lockdowns %>% filter(is.na(date)) %>% select(states) %>% as_vector()
infections %>%
    ## filter(group == "New Orleans-Metairie, LA") %>%
    filter(group %in% c("New Orleans-Metairie, LA", "Other_Louisiana")) %>%
    ggplot() + geom_line(aes(x= time, y = median), color = "deepskyblue4") + geom_ribbon(aes(x= time, ymin = lower, ymax = upper), alpha = 0.5, fill = "deepskyblue4") + facet_grid(group ~ .) + xlab("Date") + ylab("Predicted number of infections") + theme_bw()

## https://www.acpjournals.org/doi/10.7326/M20-0504
time_to_onset_gamma_shape <- 5.807
time_to_onset_gamma_rate <- 1/0.948
ndays <- 135

time_to_onset_cdf <- pgamma(seq(0, ndays, 0.1), shape=time_to_onset_gamma_shape, rate=time_to_onset_gamma_rate)
time_to_onset_pdf <- dgamma(seq(0, ndays, 0.1), shape=time_to_onset_gamma_shape, rate=time_to_onset_gamma_rate)

## Infectious period
# https://www.medrxiv.org/content/10.1101/2020.01.29.20019547v2.full.pdf
infectious_period_gamma_shape  <- 2.5
infectious_period_gamma_rate <- 0.35

infectious_period_cdf <- pgamma(seq(0, ndays, 0.1), shape=infectious_period_gamma_shape, rate=infectious_period_gamma_rate)
infectious_period_pdf <- dgamma(seq(0, ndays, 0.1), shape=infectious_period_gamma_shape, rate=infectious_period_gamma_rate)

df <- data.frame(time_to_onset = time_to_onset_pdf, infectious_period = infectious_period_pdf, n = seq(0, ndays, 0.1))

df %>%
    filter(n <= 20) %>%
    gather(key, val, -n) %>%
    ggplot(aes(n, val)) + geom_line() + facet_grid(key ~ .) + theme_bw() + ylab("Probability") + xlab("Days from infection")
ggsave("./infection_to_onset_infectious_period.svg", w= 7.5, h = 10)


## Set start at 1 and increment by 1 for each consecutive day
time_to_onset_pdf <- dgamma(1:ndays, shape=time_to_onset_gamma_shape, rate=time_to_onset_gamma_rate)
infectious_period_cdf <- pgamma(1:ndays, shape=infectious_period_gamma_shape, rate=infectious_period_gamma_rate)

## Fake data
days <- seq(1, ndays)
ncases <- c(100, rep(0, 104))
nreported <- c(rep(0, 6), rep(20, 5), rep(0, 94))
df <- tibble(
    days,
    ncases,
    nreported
)

df <- df %>%
    mutate(
        symptom_onset = df %>%
            select(days, ncases) %>%
            pmap_dbl(~{
                current_day <- .x
                df %>%
                    filter(days < current_day) %>%
                    select(days, ncases) %>%
                    pmap_dbl(~{
                        .y * (time_to_onset_pdf[current_day - .x])
                    }) %>%
                    sum()
            }),
        start_infectious = replace_na(lead(symptom_onset), 0),
        reporting_lead = lead(nreported, 5)
    )
df <- df %>%
    mutate(
        infectious = df %>%
            select(days, start_infectious) %>%
            pmap_dbl(~{
                current_day <- .x
                df %>%
                    filter(days <= current_day) %>%
                    select(days, start_infectious) %>%
                    pmap_dbl(~{
                        .y * (1 - infectious_period_cdf[(current_day + 1) - .x])
                    }) %>%
                    sum()
            }),
        travel_infectious = df %>%
            select(days, start_infectious) %>%
            pmap_dbl(~{
                current_day <- .x
                df %>%
                    filter(days <= current_day) %>%
                    select(days, start_infectious, reporting_lead) %>%
                    pmap_dbl(~{
                        travel_cases <- ifelse(..2 >= ..3, ..2 - ..3, 0) # Subtract cumulative reported cases until that day
                        travel_cases * (1 - infectious_period_cdf[(current_day + 1) - ..1])
                    }) %>%
                    sum()
            })
    )

## Plot ncases, infectious and sympton onset
time_to_onset_pdf <- dgamma(1:ndays, shape=time_to_onset_gamma_shape, rate=time_to_onset_gamma_rate)
infectious_period_cdf <- pgamma(1:ndays, shape=infectious_period_gamma_shape, rate=infectious_period_gamma_rate)

keys <- c("ncases", "nreported", "reporting_lead", "symptom_onset", "start_infectious", "infectious", "travel_infectious")
df %>%
    gather(key, val, -days) %>%
    mutate(
        key = fct_relevel(key, keys)
    ) %>%
    ggplot(aes(days, val, fill = key)) + geom_col() + facet_grid(key ~ .) + theme_bw()
ggsave("./test_symptom_onset_infectious.pdf", w=7.5, h =10)

## Infectious travellers for US states
merged_cases <- left_join(infections, cases, by=c("time" = "date", "group" = "name")) %>%
    mutate(
        confirmed_numIncrease = replace_na(confirmed_numIncrease, 0)
    )

merged_cases %>% select(time, group, median, lower, upper, confirmed)

merged_cases <- merged_cases %>%
    group_by(group) %>%
    group_modify(~{
        df <- .x
        df %>%
            mutate(
                symptom_onset = df %>%
                    select(time, median) %>%
                    pmap_dbl(~{
                        current_day <- .x
                        df %>%
                            filter(time < current_day) %>%
                            select(time, median) %>%
                            pmap_dbl(~{
                                .y * (time_to_onset_pdf[current_day - .x])
                            }) %>%
                            sum()
                    }),
                symptom_onset_lower = df %>%
                    select(time, lower) %>%
                    pmap_dbl(~{
                        current_day <- .x
                        df %>%
                            filter(time < current_day) %>%
                            select(time, lower) %>%
                            pmap_dbl(~{
                                .y * (time_to_onset_pdf[current_day - .x])
                            }) %>%
                            sum()
                    }),
                symptom_onset_upper = df %>%
                    select(time, upper) %>%
                    pmap_dbl(~{
                        current_day <- .x
                        df %>%
                            filter(time < current_day) %>%
                            select(time, upper) %>%
                            pmap_dbl(~{
                                .y * (time_to_onset_pdf[current_day - .x])
                            }) %>%
                            sum()
                    }),
                start_infectious = replace_na(lead(symptom_onset), 0),
                start_infectious_lower = replace_na(lead(symptom_onset_lower), 0),
                start_infectious_upper = replace_na(lead(symptom_onset_upper), 0),
                reporting_lead = lead(confirmed_numIncrease, 5)
            )
    }) %>%
    ungroup()

merged_cases <- merged_cases %>%
    group_by(group) %>%
    group_modify(~{
        df <- .x
        .x %>%
            mutate(
                infectious = df %>%
                    select(time, start_infectious) %>%
                    pmap_dbl(~{
                        current_day <- .x
                        df %>%
                            filter(time <= current_day) %>%
                            select(time, start_infectious) %>%
                            pmap_dbl(~{
                                .y * (1 - infectious_period_cdf[(current_day + 1) - .x])
                            }) %>%
                            sum()
                    }),
                travel_infectious = df %>%
                    select(time, start_infectious) %>%
                    pmap_dbl(~{
                        current_day <- .x
                        df %>%
                            filter(time <= current_day) %>%
                            select(time, start_infectious, reporting_lead) %>%
                            pmap_dbl(~{
                                travel_cases <- ifelse(..2 >= ..3, ..2 - ..3, 0) # Subtract cumulative reported cases until that day
                                travel_cases * (1 - infectious_period_cdf[(current_day + 1) - ..1])
                            }) %>%
                            sum()
                    }),
                infectious_lower = df %>%
                    select(time, start_infectious_lower) %>%
                    pmap_dbl(~{
                        current_day <- .x
                        df %>%
                            filter(time <= current_day) %>%
                            select(time, start_infectious_lower) %>%
                            pmap_dbl(~{
                                .y * (1 - infectious_period_cdf[(current_day + 1) - .x])
                            }) %>%
                            sum()
                    }),
                travel_infectious_lower = df %>%
                    select(time, start_infectious_lower) %>%
                    pmap_dbl(~{
                        current_day <- .x
                        df %>%
                            filter(time <= current_day) %>%
                            select(time, start_infectious_lower, reporting_lead) %>%
                            pmap_dbl(~{
                                travel_cases <- ifelse(..2 >= ..3, ..2 - ..3, 0) # Subtract cumulative reported cases until that day
                                travel_cases * (1 - infectious_period_cdf[(current_day + 1) - ..1])
                            }) %>%
                            sum()
                    }),
                infectious_upper = df %>%
                    select(time, start_infectious_upper) %>%
                    pmap_dbl(~{
                        current_day <- .x
                        df %>%
                            filter(time <= current_day) %>%
                            select(time, start_infectious_upper) %>%
                            pmap_dbl(~{
                                .y * (1 - infectious_period_cdf[(current_day + 1) - .x])
                            }) %>%
                            sum()
                    }),
                travel_infectious_upper = df %>%
                    select(time, start_infectious_upper) %>%
                    pmap_dbl(~{
                        current_day <- .x
                        df %>%
                            filter(time <= current_day) %>%
                            select(time, start_infectious_upper, reporting_lead) %>%
                            pmap_dbl(~{
                                travel_cases <- ifelse(..2 >= ..3, ..2 - ..3, 0) # Subtract cumulative reported cases until that day
                                travel_cases * (1 - infectious_period_cdf[(current_day + 1) - ..1])
                            }) %>%
                            sum()
                    })
            )
    }) %>%
    ungroup()

merged_cases %>%
    ungroup() %>%
    select(!matches("breaks|geometry|sub_parts")) %>%
    write_csv("./infectious_travellers.csv")
