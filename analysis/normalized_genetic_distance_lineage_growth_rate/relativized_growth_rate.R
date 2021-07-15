library(tidyverse)
library(data.table)
library(KFAS)
library(zoo)
library(jsonlite)
library(sf)
library(rmapshaper)
library(patchwork)
library(lubridate)
library(ggrepel)

## Read in JSON (obtained via outbreak.info from GISAID). Takes ~8 minutes
start_time <- Sys.time()
all_data <- fromJSON("./data/new_api_data.json.gz")
print(paste0(Sys.time() - start_time, " minutes"))


us_data <- all_data %>%
    as_tibble() %>%
    filter(country_id == "USA")

## Get specifi variants
## B.1.351 in South Africa
sa <- all_data %>%
    as_tibble() %>%
    filter(country == "South Africa")

## B117 in UK
uk <- all_data %>%
    as_tibble() %>%
    filter(country == "United Kingdom")

## Brazil P.1
brazil <- all_data %>%
    as_tibble() %>%
    filter(country == "Brazil")

## Massacheusetts
ma <- all_data %>%
    as_tibble() %>%
    filter(country_id == "USA" & division_id == "MA") %>%
    mutate(
        location_id = division_id,
        location = division
    )

## Cruise ship: Diamond Princess
## diamond_princess <- read_tsv("./data/gisaid_hcov-19_2021_04_05_18_diamond_princess_patient_metadata.tsv")
## diamond_princess %>%
##     group_by(`Collection date`) %>%
##     count() %>%
##     ggplot(aes(`Collection date`, n)) + geom_col() + theme_bw()

other_variants <- bind_rows(sa, uk, uk) %>%
    mutate(
        location_id = country_id,
        location = country,
        division = country,
        division_id = country_id
    ) %>%
    bind_rows(ma)

other_variants_counts <- other_variants %>%
    group_by(location, location_id, pangolin_lineage, date_collected) %>%
    count()

other_variants_counts %>%
    mutate(
        date_collected = parse_date(date_collected)
    ) %>%
    ggplot(aes(date_collected, n)) + geom_col() + facet_grid(location_id ~ pangolin_lineage, scales="free") + theme_bw()

## Read shapefile
us_county_map <- read_sf("./data/geo/tl_2019_us_county.shp") %>%
    mutate(
        location_id = paste0(STATEFP, COUNTYFP)
    ) %>%
    ms_simplify()

us_county_data <- us_data %>%
    as_tibble() %>%
    filter(location_id != "None")


us_county_counts <- us_county_data %>%
    group_by(location, location_id, division, division_id) %>%
    count()

merged_county <- us_county_map %>%
    left_join(us_county_counts, by="location_id")

## Focus on continental states
exclude_states <- c("03", "07", "14", "43", "52", "15", "72", "66", "02", "60", "69", "78")
merged_county %>%
    filter(!(STATEFP %in% exclude_states)) %>%
    ggplot() +
    geom_sf(aes(fill = n), color="#333333") +
    coord_sf(crs = "+proj=aea +lat_1=25 +lat_2=50 +lon_0=-100") +
    scale_fill_continuous(low="red", high="yellow", na.value="#ECECEC", trans="log")
    theme_void()
ggsave("./plots/us_county_sequence_sampling.pdf", w= 7.5, h =5)

## Filter down to counties with >1000 sequences
select_counties <- us_county_counts %>%
    filter(n >= 1000) %>%
    ungroup() %>%
    select(location_id) %>%
    as_vector()

## Get sequences from Louisiana
la_counties <- merged_county %>% filter(STATEFP == "22") %>% pull(location_id)
new_orleans_counties <- c("22051", "22071", "22075", "22087", "22089", "22093", "22095", "22103") # METRO_35380
shreveport_counties <- c("22015", "22017", "22031") # METRO_43340

select_counties <- c(select_counties, new_orleans_counties, shreveport_counties)

us_county_filtered <- us_county_data %>%
    filter(location_id %in% select_counties)

## Set New Orleans, Shreveport and Other Louisiana locations
us_county_filtered <- us_county_filtered %>%
    mutate(
        location = ifelse(location_id %in% new_orleans_counties, "New Orleans metro", location),
        location = ifelse(location_id %in% shreveport_counties, "Shreveport metro", location),
        location_id = ifelse(location_id %in% new_orleans_counties, "METRO_35380", location_id),
        location_id = ifelse(location_id %in% shreveport_counties, "METRO_43340", location_id)
    )

## Join other variants
us_county_filtered <- bind_rows(us_county_filtered, other_variants)

lineage_counts <- us_county_filtered %>%
    mutate(
        date_collected = parse_date(date_collected)
    ) %>%
    filter(pangolin_lineage != "None") %>%
    group_by(division, division_id, location, location_id, pangolin_lineage, date_collected) %>%
    count() %>%
    group_by(division, division_id, location_id, location, pangolin_lineage) %>%
    group_modify(~{
        .x %>%
            complete(date_collected = seq(min(date_collected)-7, max(date_collected)+7, by="day"), fill = list(n = 0))
    })

prevalence_data <- lineage_counts %>%
    group_by(division, location, division_id, location_id) %>%
    group_modify(~{
        loc_id <- .y$location_id
        ## print(loc_id)
        tmp <- .x %>%
            group_by(pangolin_lineage) %>%
            group_modify(~{ #Compute rolling lineage count
                .x %>%
                    arrange(date_collected) %>%
                    mutate(
                        lineage_count_rolling =
                            rollapply(n,FUN=function(x) mean(x, na.rm=T), width=7,align='center',fill=0)
                    )
            }) %>%
            group_by(date_collected) %>%
            group_modify(~{
                .x %>%
                    mutate(
                        prevalence = lineage_count_rolling/sum(lineage_count_rolling, na.rm=T)
                    )
            }) %>%
            group_by(pangolin_lineage) %>%
            group_modify(~{
                .x %>%
                    mutate(
                        ndays_over_threshold = length(which(prevalence >= 0.05)),
                        total_lineage_count = sum(n)
                    )
            }) %>%
            ungroup()
        if(loc_id %in% c("METRO_35380", "METRO_43340")){ # Excldue new orleans metros
            tmp <- tmp %>%
                mutate(
                    pangolin_lineage = ifelse(
                                           total_lineage_count > 50,
                                           pangolin_lineage,
                                           "Other"
                                       )
                    ) #Set lineages as "Other" if they don't match condition
        } else {
            tmp <- tmp %>%
                mutate(
                    pangolin_lineage = ifelse(
                        ndays_over_threshold >= 30,
                        pangolin_lineage,
                        "Other"
                    )
                ) #Set lineages as "Other" if they don't match condition
        }
        tmp %>%
            group_by(date_collected, pangolin_lineage) %>%
            summarise(
                n = sum(n, na.rm=T) #Recompute lineage count with new lineage "Other"
            ) %>%
            group_by(pangolin_lineage) %>%
            group_modify(~{ #Recompute rolling lineage count
                .x %>%
                    arrange(date_collected) %>%
                    mutate(
                        lineage_count_rolling = rollapply(n,FUN=function(x) mean(x, na.rm=T), width=7,align='center',fill=0)
                    )
            }) %>%
            group_by(date_collected) %>%
            group_modify(~{
                .x %>%
                    mutate(
                        prevalence = lineage_count_rolling/sum(lineage_count_rolling, na.rm=T)
                    )
            })
    })

## Check if rolling is working correctly
## prevalence_data %>%
##     ungroup() %>%
##     arrange(date_collected) %>%
##     filter(location == "New Orleans metro") %>%
##     ggplot(aes(date_collected, prevalence)) + geom_area(aes(fill = pangolin_lineage)) + facet_grid(location ~ .) + theme_bw()
    ## select(n, lineage_count_rolling, date_collected) %>%
    ## ggplot() + geom_col(aes(date_collected, n), fill = "indianred") + geom_line(aes(date_collected, lineage_count_rolling), color = "#000000") + theme_bw()

## ## Check if prevalence is working correctly
## prevalence_data %>%
##     ungroup() %>%
##     arrange(date_collected) %>%
##     filter(location == "New Orleans metro" & pangolin_lineage=="Other") %>%
##     ggplot() + geom_line(aes(date_collected, prevalence), color = "#000000") + theme_bw()

## prevalence_data %>%
##     filter(location == "United Kingdom" & pangolin_lineage == "B.1.1.7") %>%
##     ggplot(aes(date_collected, prevalence)) + geom_col() + theme_bw()

prevalence_data %>%
    ungroup() %>%
    ggplot(aes(date_collected, prevalence)) + geom_area(aes(fill = pangolin_lineage)) + facet_grid(location ~ .) + theme_bw()
ggsave("./plots/prevalence_select_counties.pdf", h = 25, w = 15)

## ny_prop <- read_tsv("./data/california_2021-03-21.tsv", col_names=c("date", "total_count", "lineage_count", "lineage", "prevalence", "prevalence_rolling"))

nbss <- function(N, seasonal = T){
# Original code form https://github.com/nicholasdavies/newcovid/blob/master/relativized_growth_rate_and_IPO_analysis_VOC_within_B117.R
    nb_ss_model <- function(N,dispersion){
        model <- SSModel(N~SSMtrend(2, Q=list(0, NA),
                                    P1=diag(c(10, 1)),
                                    a1=c(0, 0),
                                    state_names=c("abundance", "growth_rate")) +
                             SSMseasonal(7)
                        ,
                         u=rep(exp(dispersion),length(N)),distribution='negative binomial')
        fit <- fitSSM(model,inits = 0,method='L-BFGS-B',control=list(maxit=200))
    }

    if(!seasonal){
        nb_ss_model <- function(N,dispersion){
            model <- SSModel(N~SSMtrend(2, Q=list(0, NA),
                                        P1=diag(c(10, 1)),
                                        a1=c(0, 0),
                                        state_names=c("abundance", "growth_rate")),
                             u=rep(exp(dispersion),length(N)),distribution='negative binomial')
            fit <- fitSSM(model,inits = 0,method='L-BFGS-B',control=list(maxit=200))
        }
    }

    nb_log_likelihood <- function(N, disp){
        fit <- nb_ss_model(N, disp)
        ll <- logLik(fit$model, marginal = TRUE)
        return(-ll)
    }

    ### maximum marginal likelihood estimation
    mldisp <- optim(c(-1), function(disp) nb_log_likelihood(N, disp),
                    method="Brent", lower=-2, upper=2)
    fit <- nb_ss_model(N, mldisp$par)


    estimates <- KFS(fit$model, smoothing="state")
    out <- data.table(q2.5_abundance= c(qnorm(0.025, estimates$alphahat[,1], (sqrt(estimates$V[1,1,])))),
                      q97.5_abundance= c(qnorm(0.975, estimates$alphahat[,1],(sqrt(estimates$V[1,1,])))),
                      abundance = (c(estimates$alphahat[,1])),
                      q2.5_growth_rate = c(qnorm(0.025, estimates$alphahat[,2], (sqrt(estimates$V[2,2,])))),
                      q97.5_growth_rate = c(qnorm(0.975, estimates$alphahat[,2], (sqrt(estimates$V[2,2,])))),
                      q25_growth_rate = c(qnorm(0.25, estimates$alphahat[,2], (sqrt(estimates$V[2,2,])))),
                      q75_growth_rate = c(qnorm(0.75, estimates$alphahat[,2], (sqrt(estimates$V[2,2,])))),
                      growth_rate = (c(estimates$alphahat[,2])),
                      dispersion=mldisp$par[1],
                      sd_growth_rate = sqrt(estimates$V[2,2,]),
                      N=N)
    return(out)
}

## Validation
logistic_growth <- function(ndays,x0,r, N){
    round(N/(1+(((N/x0) - 1) * exp( -1 * r * ndays))))
}

cum_cases <- logistic_growth(seq(1,200), 1, 0.1, 10000)
ncases <- diff(cum_cases)

out <- nbss(ncases) %>%
    as_tibble() %>%
    mutate(
        ncases = ncases,
        norm_ncases = ncases/max(ncases),
        ndays = seq(1, length(ncases))
    )

out %>%
    ggplot() +
    geom_col(aes(ndays, norm_ncases), fill = "lightgray") +
    geom_line(aes(ndays, growth_rate)) +
    geom_ribbon(aes(ndays, ymin=q2.5_growth_rate, ymax = q97.5_growth_rate), alpha = 0.3, fill = "indianred") +
    theme_bw()
ggsave("./plots/cumulative_cases_growth_rate.png", w= 7.5, h =5)

## Cyclic trend
cyc_ncases <- 20 * sin(7 * seq(1,length(ncases)))
cyc_ncases <- round(ifelse(cyc_ncases < 0, 0, cyc_ncases))

out <- nbss(cyc_ncases) %>%
    as_tibble() %>%
    mutate(
        ncases = cyc_ncases,
        norm_ncases = cyc_ncases/max(cyc_ncases),
        ndays = seq(1, length(cyc_ncases))
    )

out <- nbss(cyc_ncases, seasonal = F) %>%
    as_tibble() %>%
    mutate(
        ncases = cyc_ncases,
        norm_ncases = cyc_ncases/max(cyc_ncases),
        ndays = seq(1, length(cyc_ncases))
    ) %>%
    inner_join(out, by="ndays", suffix = c("_seasonal", "_noseasonal"))

out %>%
    ggplot() +
    geom_col(aes(ndays, norm_ncases_seasonal), fill = "lightgray") +
    geom_line(aes(ndays, growth_rate_seasonal), color="indianred") +
    geom_ribbon(aes(ndays, ymin=q2.5_growth_rate_seasonal, ymax = q97.5_growth_rate_seasonal), alpha = 0.3, fill = "indianred") +
    geom_line(aes(ndays, growth_rate_noseasonal), color="steelblue") +
    geom_ribbon(aes(ndays, ymin=q2.5_growth_rate_noseasonal, ymax = q97.5_growth_rate_noseasonal), alpha = 0.3, fill = "steelblue") +
    theme_bw()

out <- ny_prop %>%
    filter(date <= "2021-03-06") %>%
    group_by(lineage) %>%
    group_modify(~{
        tmp <- .x %>%
            arrange(date) %>%
            mutate(
                lineage_count_rolling = rollapply(lineage_count,FUN=function(x) as.integer(mean(x, na.rm=T)), width=7,align='center',fill=NA)
            )
        nbss(tmp$lineage_count_rolling) %>%
            as_tibble() %>%
            cbind(tmp)
    })

out <- out %>%
    group_by(lineage) %>%
    group_modify(~{
        .x %>%
            mutate(
                relativized_growth_rate = (growth_rate-mean(growth_rate,na.rm=T))/sd(growth_rate,na.rm=T),
                relativized_growth_rate_q2.5 = (q2.5_growth_rate-mean(q2.5_growth_rate,na.rm=T))/sd(q2.5_growth_rate,na.rm=T),
                relativized_growth_rate_q97.5 = (q97.5_growth_rate-mean(q97.5_growth_rate,na.rm=T))/sd(q97.5_growth_rate,na.rm=T),
                lineage_count_prop = lineage_count/sum(lineage_count, na.rm=T)
            )
    })

## loc = "new_york"
loc="cali"
out %>%
    filter(lineage != "other" & date >= "2021-01-01") %>%
    ggplot(aes(date, relativized_growth_rate, color = lineage)) +
    geom_line() +
    theme_bw()
ggsave(paste0("./plots/",loc,"_growth_rate.png"), w= 7.5, h =5)

out %>%
    ggplot(aes(date, lineage_count, fill = lineage)) +
    geom_col() +
    geom_line(aes(date, lineage_count_rolling)) +
    theme_bw() +
    facet_grid(lineage ~ ., scales="free")
ggsave(paste0("./plots/", loc,"_lineage_counts.png"), w= 5, h =15)

out %>%
    filter(date >= "2021-02-01") %>%
    ggplot(aes(lineage, relativized_growth_rate)) + geom_violin() + geom_boxplot(width=0.1) + theme_bw()
ggsave(paste0("./plots/", loc,"_last_month.png"), w= 5, h =15)

## Get growth rates of New Orleans
lineage_growth <- prevalence_data %>%
    filter(pangolin_lineage != "Other") %>%
    group_by(division, division_id, location_id, location, pangolin_lineage) %>%
    group_modify(~{
        ## First week with at least 5 cases
        tmp <- .x %>%
            arrange(date_collected) %>%
            mutate(
                date_collected_epi_week = epiweek(date_collected),
                date_collected_year = year(date_collected)
            )

        first_week <- tmp %>%
            group_by(date_collected_year, date_collected_epi_week) %>%
            summarise(
                nday_threshold = length(which(n > 0))
            ) %>%
            filter(nday_threshold >= 5) %>%
            arrange(date_collected_year, date_collected_epi_week) %>%
            head(1)
        res <- data.frame()

        if(nrow(first_week) > 0){

            last_week <- tmp %>%
                filter(
                    date_collected_epi_week >= first_week %>% pull(date_collected_epi_week) &
                    date_collected_year >= first_week %>% pull(date_collected_year)
                ) %>%
                group_by(date_collected_year, date_collected_epi_week) %>%
                summarise(
                    nday_threshold = length(which(n == 0))
                ) %>%
                filter(nday_threshold >= 5) %>%
                arrange(date_collected_year, date_collected_epi_week) %>%
                head(1)

            start_date <- tmp %>%
                filter(
                    date_collected_epi_week == first_week %>% pull(date_collected_epi_week) &
                    date_collected_year == first_week %>% pull(date_collected_year) &
                    n > 0
                ) %>%
                arrange(date_collected) %>%
                pull(date_collected) %>%
                first()
            end_date <- NA
            if(nrow(last_week) > 0){
                end_date <- tmp %>%
                    filter(
                        date_collected_epi_week == last_week %>% pull(date_collected_epi_week) &
                        date_collected_year == last_week %>% pull(date_collected_year) &
                        n == 0
                    ) %>%
                    arrange(date_collected) %>%
                    pull(date_collected) %>%
                    first()
            } else {
                end_date = tmp %>% arrange(date_collected) %>% pull(date_collected) %>% last()
            }

            ndays <- end_date - start_date

            print(paste(.y$location, .y$pangolin_lineage, start_date, end_date, collapse=" "))

            tmp <- .x %>%
                filter(date_collected >= start_date & date_collected <= end_date) %>%
                mutate(
                    nday = date_collected - start_date
                )

            if(tmp %>% pull(n) %>% sum(na.rm=TRUE) >= 10){ #At leas 10 sequences across whole time period
                growth_df <- tmp %>%
                    pull(lineage_count_rolling) %>%
                    round() %>%
                    nbss() %>%
                    as_tibble() %>%
                    bind_cols(tmp)

                p1 <- growth_df %>%
                    ggplot() +
                    geom_col(aes(date_collected, n), fill="indianred") +
                    geom_line(aes(date_collected, lineage_count_rolling), color="#000000") +
                    scale_x_date(limits=c(start_date-3, start_date + ndays+3)) +
                    theme_bw()

                p2 <- growth_df %>%
                    ggplot() +
                    geom_line(aes(date_collected, growth_rate), color="indianred") +
                    geom_ribbon(aes(date_collected, ymin = q2.5_growth_rate, ymax = q97.5_growth_rate), color="#000000", alpha = 0.2) +
                    scale_x_date(limits=c(start_date-3, start_date + ndays + 3)) +
                    scale_y_continuous(limits=c(
                                           growth_df %>% select(q2.5_growth_rate) %>% pull %>% min - 0.5,
                                           growth_df %>% select(q97.5_growth_rate) %>% pull %>% max + 0.5
                                       )) +
                    theme_bw()

                p1/p2
                ggsave(paste0("./plots/growth_rates/", .y$division, "_", .y$location, "_", .y$pangolin_lineage,"_growth.pdf"), w = 7.5, h =10)

                res <- growth_df
            }
        }
        res
    })

lineage_growth %>%
    write_csv("./lineage_growth.csv")

mean_lineage_growth <- lineage_growth %>%
    mutate(
        location_lineage = paste0(division, "_", location, "_", pangolin_lineage)
    ) %>%
    group_by(location_lineage, location, division, location_id, division_id, pangolin_lineage) %>%
    group_modify(~{
        end_date <- .x %>%
            arrange(nday) %>%
            filter(lineage_count_rolling == max(lineage_count_rolling)) %>%
            pull(nday) %>%
            first()
        if(length(end_date) == 0){
            end_date <- .x %>% arrange(nday) %>% pull(nday) %>% max()
        }
        end_date <- 10
        .x %>%
            filter(nday <= end_date) %>%
            summarise(mean_growth_rate = mean(growth_rate), nday_threshold = max(nday), first_date = min(date_collected), end_date = min(date_collected) +end_date)
    })

tmp <- inner_join(mean_lineage_growth, us_county_filtered, by=c("location_id", "division_id", "pangolin_lineage", "location", "division"))

## Calculate prevalence
prev_10_day <- tmp %>%
    filter(date_collected >= first_date & date_collected <= end_date) %>%
    group_by(location, division, location_id, division_id) %>%
    group_modify(~{
        tmp <- .x
        tmp %>%
            group_by(pangolin_lineage) %>%
            group_modify(~{
                .x %>%
                    summarise(
                        prev_10_day = n()/nrow(tmp)
                    )
            })
    }) %>%
    arrange(-prev_10_day)

## Merge
inner_join(prev_10_day, mean_lineage_growth, by=c("location_id", "division_id", "pangolin_lineage", "location", "division")) %>%
    write_csv("./mean_lineage_growth.csv")

## Get accession ids
tmp <- inner_join(mean_lineage_growth, us_county_filtered, by=c("location_id", "division_id", "pangolin_lineage", "location", "division"))

tmp %>%
    filter(date_collected >= first_date & date_collected <= end_date) %>%
    select(accession_id, location, division, division_id, first_date, end_date, date_collected) %>%
    write_csv("./accessin_ids.csv")


#Check if filter worked
tmp %>%
    group_by(location.x, division.x, pangolin_lineage) %>%
    summarise(start_date = first(first_date), min_date = min(date_collected), max_date = max(date_collected), end_date = first(end_date)) %>%
    ungroup() %>%
    filter(
        start_date < min_date | max_date > end_date
    ) #Should be 0 rows


## Plot mean growth rate vs IPO
mean_lineage_growth %>%
    filter(nday_threshold == 10) %>%
    mutate(
        label_point = ifelse(mean_growth_rate > 0.08 | str_detect(location_lineage, "United Kingdom|South Africa"), location_lineage, NA)
    ) %>%
    ggplot(aes(first_date, mean_growth_rate)) +
    geom_point() +
    geom_text_repel(aes(first_date, mean_growth_rate, label = label_point)) +
    theme_bw()
ggsave("./plots/mean_growth_rate_vs_first_date.pdf", w=15, h =7.5)

lineage_growth %>%
    filter(pangolin_lineage == "B.1") %>%
    mutate(
        location_lineage = paste0(location, pangolin_lineage)
    ) %>%
    ggplot(aes(nday, growth_rate, color = location_lineage)) + geom_line() + theme_bw()
ggsave("./plots/growth_B.1.pdf", w= 7.5, h =5)

lineage_growth %>%
    filter(location %in% c("United Kingdom", "South Africa", "New Orleans metro") & pangolin_lineage %in% c("B.1", "B.1.1.7", "B.1.351", "B.1.177")) %>%
    mutate(
        location_lineage = paste0(location, pangolin_lineage)
    ) %>%
    ggplot(aes(nday, growth_rate, color = location_lineage)) + geom_line() + theme_bw()
ggsave("./plots/growth_voc.pdf", w= 7.5, h =5)


lineage_growth %>%
    filter(location %in% c("New Orleans metro", "New York", "Massachusetts") & pangolin_lineage %in% c("B.1")) %>%
    mutate(
        location_lineage = paste0(location, pangolin_lineage)
    ) %>%
    ggplot(aes(nday, growth_rate, color = location_lineage)) + geom_line() + theme_bw()
ggsave("./plots/growth_ss.pdf", w= 7.5, h =5)

dominant_lineages <- lineage_growth %>%
    group_by(pangolin_lineage, location) %>%
    summarise(nday_threshold = length(which(prevalence > 0.5))) %>%
    filter(nday_threshold >= 10)

lineage_growth %>%
    inner_join(dominant_lineages, by=c("pangolin_lineage", "location")) %>%
    filter(!(location %in% c("South Africa", "United Kingdom"))) %>%
    mutate(
        location_lineage = paste0(location, pangolin_lineage)
    ) %>%
    ggplot(aes(nday, growth_rate, color = location_lineage)) + geom_line() + theme_bw()
ggsave("./plots/growth_dominant.pdf", w= 7.5, h =5)

nola_combined <- lineage_growth %>%
    filter(location == "New Orleans metro") %>%
    group_by(date_collected) %>%
    summarise(lineage_count_rolling = round(sum(lineage_count_rolling, na.rm=T)))

nola_growth <- nbss(nola_combined$lineage_count_rolling) %>%
    as_tibble() %>%
    cbind(nola_combined)

nola_growth %>%
    ggplot(aes(date_collected, growth_rate)) +
    geom_line(aes(date_collected, growth_rate), color="indianred") +
    geom_ribbon(aes(date_collected, ymin = q2.5_growth_rate, ymax = q97.5_growth_rate), color="#000000", alpha = 0.2)

nola_growth <- nola_growth %>%
    mutate(
        nday = date_collected - min(date_collected),
        location = "New Orleans metro",
        pangolin_lineage = "B.1 + B.1.595",
        prevalence = 1
    )

tmp <- lineage_growth %>%
    filter(location != "New Orleans metro") %>%
    bind_rows(nola_growth)

## Plot B.1 and B.1.595 together
tmp %>%
    filter(location %in% c("United Kingdom", "South Africa", "New Orleans metro") & pangolin_lineage %in% c("B.1", "B.1.1.7", "B.1.351", "B.1.177", "B.1 + B.1.595")) %>%
    mutate(
        location_lineage = paste0(location, pangolin_lineage)
    ) %>%
    ggplot(aes(nday, growth_rate, color = location_lineage)) + geom_line() + theme_bw()
ggsave("./plots/growth_voc_nola_combined.pdf", w= 7.5, h =5)


tmp %>%
    filter(location %in% c("New Orleans metro", "New York", "Massachusetts") & pangolin_lineage %in% c("B.1", "B.1 + B.1.595")) %>%
    mutate(
        location_lineage = paste0(location, pangolin_lineage)
    ) %>%
    ggplot(aes(nday, growth_rate, color = location_lineage)) + geom_line() + theme_bw()
ggsave("./plots/growth_ss_nola_combined.pdf", w= 7.5, h =5)

dominant_lineages <- tmp %>%
    group_by(pangolin_lineage, location) %>%
    summarise(nday_threshold = length(which(prevalence > 0.5))) %>%
    filter(nday_threshold >= 10)

tmp %>%
    inner_join(dominant_lineages, by=c("pangolin_lineage", "location")) %>%
    filter(!(location %in% c("South Africa", "United Kingdom"))) %>%
    mutate(
        location_lineage = paste0(location, pangolin_lineage)
    ) %>%
    ggplot(aes(nday, growth_rate, color = location_lineage)) + geom_line() + theme_bw()
ggsave("./plots/growth_dominant_nola_combined.pdf", w= 7.5, h =5)
