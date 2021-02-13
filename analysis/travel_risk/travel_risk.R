library(tidyverse)
library(lubridate)
library(rmapshaper)
library(sf)
library(outbreakinfo)
library(gganimate)

merged_cases <- read_csv("../estimated_infections/infectious_travellers.csv", col_types = cols(
  group = col_character(),
  time = col_datetime(format = ""),
  median = col_double(),
  lower = col_double(),
  upper = col_double(),
  `_id` = col_character(),
  `_score` = col_double(),
  admin_level = col_double(),
  confirmed = col_double(),
  confirmed_doublingRate = col_double(),
  confirmed_firstDate = col_date(format = ""),
  confirmed_newToday = col_logical(),
  confirmed_numIncrease = col_double(),
  confirmed_numIncrease_per_100k = col_double(),
  confirmed_pctIncrease = col_double(),
  confirmed_per_100k = col_double(),
  confirmed_rolling = col_double(),
  confirmed_rolling_per_100k = col_double(),
  country_gdp_per_capita = col_double(),
  country_iso3 = col_character(),
  country_name = col_character(),
  country_population = col_double(),
  dead = col_double(),
  dead_firstDate = col_date(format = ""),
  dead_newToday = col_logical(),
  dead_numIncrease = col_double(),
  dead_numIncrease_per_100k = col_double(),
  dead_per_100k = col_double(),
  dead_rolling = col_double(),
  dead_rolling_per_100k = col_double(),
  `first_dead-first_confirmed` = col_double(),
  gdp_last_updated = col_double(),
  iso3 = col_character(),
  lat = col_double(),
  location_id = col_character(),
  long = col_double(),
  mostRecent = col_logical(),
  population = col_double(),
  recovered = col_double(),
  recovered_firstDate = col_logical(),
  recovered_newToday = col_logical(),
  recovered_numIncrease = col_double(),
  recovered_numIncrease_per_100k = col_double(),
  recovered_per_100k = col_double(),
  recovered_rolling = col_double(),
  recovered_rolling_per_100k = col_double(),
  wb_region = col_character(),
  confirmed_rolling_14days_ago = col_double(),
  confirmed_rolling_14days_ago_diff = col_double(),
  confirmed_rolling_14days_ago_diff_per_100k = col_double(),
  confirmed_rolling_14days_ago_per_100k = col_double(),
  dead_rolling_14days_ago = col_double(),
  dead_rolling_14days_ago_diff = col_double(),
  dead_rolling_14days_ago_diff_per_100k = col_double(),
  dead_rolling_14days_ago_per_100k = col_double(),
  recovered_rolling_14days_ago = col_double(),
  recovered_rolling_14days_ago_diff = col_double(),
  recovered_rolling_14days_ago_diff_per_100k = col_double(),
  recovered_rolling_14days_ago_per_100k = col_double(),
  daysSince100Cases = col_double(),
  dead_doublingRate = col_double(),
  dead_pctIncrease = col_double(),
  testing_checkTimeEt = col_datetime(format = ""),
  testing_commercialScore = col_double(),
  testing_dataQualityGrade = col_character(),
  testing_dateChecked = col_datetime(format = ""),
  testing_dateModified = col_datetime(format = ""),
  testing_death = col_double(),
  testing_deathIncrease = col_double(),
  testing_fips = col_character(),
  testing_grade = col_character(),
  testing_hospitalizedIncrease = col_double(),
  testing_lastUpdateEt = col_datetime(format = ""),
  testing_negative = col_double(),
  testing_negativeIncrease = col_double(),
  testing_negativeRegularScore = col_double(),
  testing_negativeScore = col_double(),
  testing_posNeg = col_double(),
  testing_positive = col_double(),
  testing_positiveIncrease = col_double(),
  testing_positiveScore = col_double(),
  testing_score = col_double(),
  testing_total = col_double(),
  testing_totalTestResults = col_double(),
  testing_totalTestResultsIncrease = col_double(),
  testing_totalTestResultsSource = col_character(),
  testing_totalTestsViral = col_double(),
  daysSince10Deaths = col_double(),
  daysSince50Deaths = col_double(),
  testing_hospitalizedCurrently = col_double(),
  testing_positiveCasesViral = col_double(),
  testing_hospitalized = col_double(),
  testing_hospitalizedCumulative = col_double(),
  testing_probableCases = col_double(),
  testing_recovered = col_double(),
  testing_deathConfirmed = col_double(),
  testing_deathProbable = col_character(),
  testing_inIcuCumulative = col_double(),
  testing_totalTestsPeopleAntibody = col_character(),
  testing_onVentilatorCumulative = col_double(),
  cbsa = col_character(),
  symptom_onset = col_double(),
  symptom_onset_lower = col_double(),
  symptom_onset_upper = col_double(),
  start_infectious = col_double(),
  start_infectious_lower = col_double(),
  start_infectious_upper = col_double(),
  reporting_lead = col_double(),
  infectious = col_double(),
  travel_infectious = col_double(),
  infectious_lower = col_double(),
  travel_infectious_lower = col_double(),
  infectious_upper = col_double(),
  travel_infectious_upper = col_double()
))

merged_cases <- merged_cases %>%
    mutate(
        time = ymd(time),
        group = str_replace(group, "_", " ")
    )

## Add population for all states in NA
merged_cases <- merged_cases %>%
    group_by(group) %>%
    group_modify(~{
        .x %>%
            arrange(time) %>%
            mutate(
                population = last(population),
                iso3 = last(iso3)
            )
    })

## exclude Alaska, Oklahoma and Montana
merged_cases <- merged_cases %>%
    filter(!(group %in% c("Alaska", "Oklahoma", "Montana")))

## IFR estimate
## Deadths with 14 day delay
cumulative_death <- getLocationData("New Orleans-Metairie, LA") %>%
    filter(date == "2020-05-29") %>%
    select(dead)

merged_cases %>%
    ungroup() %>%
    filter(group == "New Orleans-Metairie, LA") %>%
    summarise(
        ifr_median = cumulative_death/sum(median, na.rm=TRUE) * 100,
        ifr_lower = cumulative_death/sum(lower, na.rm=TRUE) * 100,
        ifr_upper = cumulative_death/sum(upper, na.rm=TRUE) * 100
    )

## Get air travel
air_states <- read_csv("./data/air_travel_sub_louisiana.csv") %>%
    mutate(
        State_name_from = str_replace(State_name_from, "_", " ")
    ) %>%
    filter(State_name_to == "New_Orleans")

air_travel_risk <- merged_cases %>%
    mutate(
        MONTH = paste0(month(time, label = TRUE, abbr=FALSE), "-2020")
    ) %>%
    inner_join(air_states, by=c("group" = "State_name_from", "MONTH" = "Date")) %>%
    mutate(
        infectious_travellers = (travel_infectious/population),
        import_risk = infectious_travellers * n,
        infectious_travellers_lower = (travel_infectious_lower/population),
        import_risk = infectious_travellers_lower * n,
        infectious_travellers_upper = (travel_infectious_upper/population),
        import_risk = infectious_travellers_upper * n
    )

air_travel_risk %>%
    select(group, time, confirmed, confirmed_rolling, population, travel_infectious, travel_infectious_lower, travel_infectious_upper, infectious_travellers, infectious_travellers_lower, infectious_travellers_upper, import_risk, n) %>%
    write_csv("./air_travel_risk.csv")

air_travel_risk %>%
    filter(time < as.Date("2020-03-15")) %>%
    select(group, time, infectious_travellers, import_risk) %>%
    gather(key, val, -group, -time) %>%
    ggplot(aes(time, val, fill = group)) + geom_col() + facet_grid(group ~ key)
ggsave("./travel_risk_till_Mar.pdf", w = 20, h = 10)

air_travel_risk %>%
    filter(time < as.Date("2020-03-01")) %>%
    select(group, time, infectious_travellers, import_risk) %>%
    gather(key, val, -group, -time) %>%
    ggplot(aes(time, val, fill = group)) + geom_col() + facet_grid(group ~ key) + geom_vline(xintercept = as.Date("2020-02-25"), color="red")
ggsave("./travel_risk_per_capita_till_Feb.pdf", w = 20, h = 10)

## Get mobility data
mobility_states <- read_csv("./data/travel_weekly_sub_louisiana.csv") %>%
    mutate(
        State_name_from = str_replace(State_name_from, "_", " ")
    ) %>%
    filter(State_name_to == "New_Orleans")

mobility_travel_risk <- merged_cases %>%
    mutate(
        WEEK = isoweek(time)
    ) %>%
    inner_join(mobility_states, by=c("group" = "State_name_from", "WEEK" = "iso_week")) %>%
    mutate(
        infectious_travellers = (travel_infectious/population),
        import_risk = infectious_travellers * n,
        infectious_travellers_lower = (travel_infectious_lower/population),
        import_risk_lower = infectious_travellers_lower * n,
        infectious_travellers_upper = (travel_infectious_upper/population),
        import_risk_upper = infectious_travellers_upper * n
    )

mobility_travel_risk %>%
    select(group, time, confirmed, confirmed_rolling, population, travel_infectious, travel_infectious_lower, travel_infectious_upper, infectious_travellers, infectious_travellers_lower, infectious_travellers_upper, import_risk, import_risk_lower, import_risk_upper, n) %>%
    write_csv("./mobility_travel_risk.csv")


mobility_travel_risk %>%
    filter(time < as.Date("2020-03-15") & group != "New Orleans") %>%
    select(group, time, infectious_travellers, import_risk) %>%
    gather(key, val, -group, -time) %>%
    ggplot(aes(time, val, fill = group)) + geom_col() + facet_grid(group ~ key)
ggsave("../mobility_travel_risk_till_Mar.pdf", w = 20, h = 10)

mobility_travel_risk %>%
    filter(time < as.Date("2020-03-01") & group != "New Orleans") %>%
    select(group, time, infectious_travellers, import_risk) %>%
    gather(key, val, -group, -time) %>%
    ggplot(aes(time, val, fill = group)) + geom_col() + facet_grid(group ~ key) + geom_vline(xintercept = as.Date("2020-02-25"), color="red")
ggsave("../mobility_travel_risk_till_Feb.pdf", w = 20, h = 10)

## Travel into every state
us_state_mobility <- read_csv("./data/us_state_travel_matrix.csv") %>%
    filter(State_name_to != State_name_from) %>%
    filter(!(State_name_from %in% c("Alaska", "Oklahoma", "Montana")) & !(State_name_to %in% c("Alaska", "Oklahoma", "Montana")))

state_code_mapping <- merged_cases %>%
    distinct(group, iso3) %>%
    ungroup() %>%
    mutate(
        iso3 = ifelse(group == "New Orleans-Metairie, LA", "New_Orleans", ifelse(group == "Shreveport-Bossier City, LA", "Shreveport_Bossier", ifelse(group == "Other Louisiana", "Other_Louisiana", iso3))),
        group = ifelse(group == "New Orleans-Metairie, LA", "New Orleans", ifelse(group == "Shreveport-Bossier City, LA", "Shreveport Bossier", group))
    )

## Add iso3 for New Orleans and Shreveport Bossier
merged_cases <- merged_cases %>%
    ungroup() %>%
    mutate(
        iso3 = ifelse(group == "New Orleans-Metairie, LA", "New_Orleans", ifelse(group == "Shreveport-Bossier City, LA", "Shreveport_Bossier", ifelse(group == "Other Louisiana", "Other_Louisiana", iso3))),
        group = ifelse(group == "New Orleans-Metairie, LA", "New Orleans", ifelse(group == "Shreveport-Bossier City, LA", "Shreveport Bossier", group))
    )

us_state_mobility <- us_state_mobility %>%
    drop_na(State_name_from) %>%
    mutate(
        State_name_from = State_name_from %>%
            map_chr(~{
                if(nchar(.x) == 2){
                    paste0("US-", .x)
                } else {
                    state_code_mapping %>%
                        filter(group == str_replace(.x, "_", " ")) %>%
                        select(iso3) %>%
                        first()
                }
            }),
        State_name_to = State_name_to %>%
            map_chr(~{
                if(nchar(.x) == 2){
                    paste0("US-", .x)
                } else {
                    state_code_mapping %>%
                        filter(group == str_replace(.x, "_", " ")) %>%
                        select(iso3) %>%
                        first()
                }
            })
    )

## Plot heatmaps for risk into each US state
us_state_mobility %>%
    distinct(State_name_from) %>%
    as_vector() %>%
    walk(~{
        tmp <- us_state_mobility %>%
            filter(State_name_to == .x)
        tmp_risk <- merged_cases %>%
            ungroup() %>%
            mutate(
                WEEK = isoweek(time)
            ) %>%
            group_by(iso3, group, WEEK) %>%
            summarise(
                travel_infectious = sum(travel_infectious),
                population = last(population)
            ) %>%
            inner_join(tmp, by=c("iso3" = "State_name_from", "WEEK" = "iso_week")) %>%
            mutate(
                infectious_travellers_per_100k = (travel_infectious/population),
                import_risk = infectious_travellers_per_100k * n
            )
        write_csv(tmp_risk, paste0("../risk_csv/", .x,".csv"))
        tmp_risk %>%
            drop_na(import_risk) %>%
            ggplot(aes(WEEK, group, fill = import_risk)) + geom_tile() + geom_text(aes(WEEK, group, label = round(import_risk, 2))) + scale_fill_distiller(palette = "YlOrRd", direction = 1) + theme_bw()
        ggsave(paste0("./plots/travel_risk_into_", .x,".pdf"), w = 20, h = 10)
        print(.x)
    })

## Plot map for export risk from NOLA. Export risk into all other states as % of risk coming into ech individual state
## us_state_mobility <- read_csv("../../export_risk/us_state_travel_matrix.csv")
## Download shapefile from https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
us_admin1 <- read_sf("./data/geo/cb_2018_us_state_5m/cb_2018_us_state_5m.shp")

export_risk_nola <- us_state_mobility %>%
    group_by(State_name_to) %>%
    group_modify(~{
        tmp_risk <- merged_cases %>%
            ungroup() %>%
            mutate(
                WEEK = isoweek(time)
            ) %>%
            group_by(iso3, group, WEEK) %>%
            summarise(
                travel_infectious = sum(travel_infectious),
                population = last(population)
            ) %>%
            inner_join(.x, by=c("iso3" = "State_name_from", "WEEK" = "iso_week")) %>%
            mutate(
                infectious_travellers_per_100k = (travel_infectious/population),
                import_risk = infectious_travellers_per_100k * n
            )
        tmp_risk %>%
            group_by(WEEK) %>%
            group_modify(~{
                .x %>%
                    mutate(
                        import_risk = (import_risk/sum(import_risk, na.rm=T)) * 100
                    ) %>%
                    filter(iso3 == "New_Orleans")
            })
    })

crs <- st_crs("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

## Filter to contiguous USA
tmp_admin1 <- us_admin1 %>%
    filter(!(NAME %in% c("Hawaii", "Puerto Rico", "Alaska", "Commonwealth of the Northern Mariana Islands", "American Samoa", "Guam", "United States Virgin Islands")))

## Animate import risk
export_risk_nola_fill_na <- export_risk_nola

statefips_list <- export_risk_nola_fill_na %>% ungroup() %>% distinct(StateFIPS_cbg) %>% mutate(import_risk = NA)

## Fill in WEEK for states with NA
export_risk_nola_fill_na <- export_risk_nola_fill_na %>%
    group_by(WEEK) %>%
    group_modify(~{
        tmp <- anti_join(statefips_list, .x, by="StateFIPS_cbg")
        bind_rows(.x, tmp)
    })

merged_map <- left_join(tmp_admin1, export_risk_nola_fill_na, by=c("STATEFP" = "StateFIPS_cbg"))

export_risk_nola_anim <- merged_map %>%
    filter(WEEK < 20) %>%
    mutate(
        import_risk = ifelse(NAME == "Louisiana", NA, import_risk)
    ) %>%
    ggplot() +
    geom_sf(aes(geometry = geometry, fill = import_risk), color="#333333", size = 0.1) +
    coord_sf(crs = crs) +
    scale_fill_distiller(palette="YlOrRd", direction = 1, limits = c(0, 60), name="Import Risk from NOLA") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) +
    transition_states(
        WEEK,
        transition_length = 2,
        state_length = 1
    ) +
    labs(title = 'Week {closest_state}') +
    enter_fade() +
    exit_fade() +
    ease_aes('linear')

animate(export_risk_nola_anim, renderer = ffmpeg_renderer(), width = 2000, height = 1400, res = 300)
anim_save("./export_risk_nola.mp4")

## Plot static images
export_risk_nola %>%
    group_by(WEEK) %>%
    group_walk(~{
        .x <- .x %>%
            filter(State_name_to != "New_Orleans" & State_name_to != "Shreveport_Bossier" & State_name_to != "Other_Louisiana") %>%
            ungroup()
        merged_map <- left_join(tmp_admin1, .x, by=c("STATEFP" = "StateFIPS_cbg"))
        ggplot() +
            geom_sf(aes(geometry = geometry, fill = import_risk), color="#333333", size = 0, data = merged_map) +
            coord_sf(crs = crs) +
            scale_fill_distiller(palette="YlOrRd", direction = 1, limits = c(0, 60)) +
            ggtitle(paste0("Week ", .y)) +
            theme_void() +
            theme(plot.title = element_text(hjust = 0.5))
        ggsave(paste0("../plot_maps/",.y,".pdf"), w = 20, h = 10)
    })

## Plot within Louisiana
## Download shapefile from https://catalog.data.gov/dataset/tiger-line-shapefile-2019-nation-u-s-current-county-and-equivalent-national-shapefile
us_admin_2 <- read_sf("./data/geo/tl_2019_us_county.shp")
la_admin_2 <- us_admin_2 %>%
    filter(STATEFP == "22")

## Shreveport and New Orleans counties
shreveport_bossier_counties <- c("015", "017", "031")
nola_counties <- c("051", "071", "075", "087", "089","093", "095", "103")

shreveport_bossier_sf <- la_admin_2 %>%
    filter(COUNTYFP %in% shreveport_bossier_counties) %>%
    st_union()

nola_sf <- la_admin_2 %>%
    filter(COUNTYFP %in% nola_counties) %>%
    st_union()

other_la_sf <- la_admin_2 %>%
    filter(!(COUNTYFP %in% c(nola_counties, shreveport_bossier_counties))) %>%
    st_union()

tmp_la_sf <- data.frame(name = c("Shreveport_Bossier", "New_Orleans", "Other_Louisiana"))
st_geometry(tmp_la_sf) <- st_sfc(st_multipolygon(shreveport_bossier_sf), st_multipolygon(nola_sf), st_multipolygon(other_la_sf))
st_crs(tmp_la_sf) <- st_crs(la_admin_2)

## Animate
la_region_list <- data.frame(State_name_to = c("Shreveport_Bossier", "Other_Louisiana", "New_Orleans")) %>%mutate(import_risk = NA)

## Fill in WEEK for regions with NA
export_risk_nola_within_la_fill_na <- export_risk_nola %>%
    filter(State_name_to %in% c("Shreveport_Bossier", "Other_Louisiana")) %>%
    group_by(WEEK) %>%
    group_modify(~{
        tmp <- anti_join(la_region_list, .x, by="State_name_to")
        bind_rows(.x, tmp)
    })

merged_map <- left_join(tmp_la_sf, export_risk_nola_within_la_fill_na, by=c("name" = "State_name_to"))

export_risk_nola_anim <- merged_map %>%
    filter(WEEK < 20) %>%
    mutate(
        import_risk = ifelse(name == "New_Orleans", NA, import_risk)
    ) %>%
    ggplot() +
    geom_sf(aes(geometry = geometry, fill = import_risk), color="#333333", size = 0.1) +
    coord_sf(crs = crs) +
    scale_fill_distiller(palette="YlOrRd", direction = 1, limits = c(0, 85), name="Import Risk from NOLA") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) +
    transition_states(
        WEEK,
        transition_length = 2,
        state_length = 1
    ) +
    labs(title = 'Week {closest_state}') +
    enter_fade() +
    exit_fade() +
    ease_aes('linear')

animate(export_risk_nola_anim, renderer = ffmpeg_renderer(), width = 2000, height = 1400, res = 300)
anim_save("./plots/export_risk_nola_within_la.mp4")


## Static images
export_risk_nola %>%
    filter(State_name_to %in% c("Shreveport_Bossier", "Other_Louisiana")) %>%
    group_by(WEEK) %>%
    group_walk(~{
        merged_map <- left_join(tmp_la_sf, .x, by=c("name" = "State_name_to"))
        ggplot() +
            geom_sf(aes(geometry = geometry, fill = import_risk), color="#333333", size = 0, data = merged_map) +
            coord_sf(crs = crs) +
            scale_fill_distiller(palette="YlOrRd", direction = 1, limits = c(0, 85)) +
            ## ggtitle(paste0("Week ", .y)) +
            theme_void() +
            theme(plot.title = element_text(hjust = 0.5))
        ggsave(paste0("../plot_maps/inset_louisiana_",.y,".pdf"), w = 10, h = 10)
    })

## Export from Shreveport Bossier
## from_nola_travel_risk <- us_state_mobility %>%
##     filter(State_name_from == "New_Orleans" & State_name_to != "New_Orleans")

from_nola_travel_risk <- us_state_mobility %>%
    filter(State_name_from == "Shreveport_Bossier" & State_name_to != "Shreveport_Bossier")

tmp_risk <- merged_cases %>%
    ungroup() %>%
    mutate(
        WEEK = isoweek(time),
        iso3 = ifelse(group == "New Orleans-Metairie, LA", "New_Orleans", ifelse(group == "Shreveport-Bossier City, LA", "Shreveport_Bossier", iso3))
    ) %>%
    group_by(iso3, group, WEEK) %>%
    summarise(
        travel_infectious = sum(travel_infectious),
        population = last(population)
    ) %>%
    inner_join(from_nola_travel_risk, by=c("iso3" = "State_name_from", "WEEK" = "iso_week")) %>%
    mutate(
        infectious_travellers_per_100k = (travel_infectious/population),
        import_risk = infectious_travellers_per_100k * n,
        WEEK = as.factor(WEEK),
        group = as.factor(group)
    )

tmp_risk %>%
    mutate(
        import_risk = replace_na(import_risk, 0)
    ) %>%
    filter(WEEK != 20) %>%
    ggplot(aes(WEEK, reorder(State_name_to, import_risk), fill = import_risk)) + geom_tile() + scale_fill_distiller(palette = "Blues", direction = 1, trans="pseudo_log", breaks = c(0, 100, 200, 300, 400, 500)) + theme_bw()
    ## ggplot(aes(WEEK, reorder(State_name_to, import_risk), fill = import_risk)) + geom_tile() + scale_fill_distiller(palette = "Blues", direction = 1, trans="pseudo_log", breaks = c(0, 500, 1000, 2000, 3000)) + theme_bw()
## ggsave("./plots/travel_risk_from_nola.svg", w= 15, h = 10)
ggsave("./plots/travel_risk_from_shreveport_bossier.svg", w= 15, h = 10)

## Air travel outgoing from NOLA
air_travel_1 <- read_csv("./data/LA_to_USA/PAX_Louisiana_to_USA_Feb-March2020.csv", col_types = cols(.default = "c"))

air_travel_2 <- read_csv("./data/LA_to_USA/PAX_Louisiana-to-USA_April-May2020.csv", col_types = cols(.default = "c")) %>%
    mutate(
        `Destination Lat` = as.numeric(`Destination Lat`),
        `Destination Long` = as.numeric(`Destination Long`)
        )

## Geo join air travel 2 to get destination province
air_travel_2 <- st_as_sf(air_travel_2, coords=c("Destination Long", "Destination Lat"), crs = st_crs(us_admin1))
air_travel_2 <- st_join(air_travel_2, us_admin1, join = st_within) %>%
    mutate(`Destination Province` = NAME) %>%
    select(Date, `Origin Airport`, `Origin City`, `Origin Country`, `Destination Airport`, `Destination City`, `Destination Country`, `Total Volume`, `Destination Province`) %>%
    as_tibble()

## Combine both air travel documents
air_travel <- bind_rows(air_travel_1, air_travel_2) %>%
    mutate(
        `Total Volume` = str_replace(`Total Volume`, ",", "") %>% as.integer()
    )

air_travel_from_nola <- air_travel %>%
    filter(`Origin City` == "New Orleans") %>%
    group_by(`Destination Province`, Date) %>%
    summarise(n = sum(`Total Volume`)) %>%
    drop_na(n)

nola_cases <- merged_cases %>%
    filter(iso3 == "New_Orleans") %>%
    mutate(
        MONTH = paste0(month(time, label=T, abbr=F), "-2020")
    ) %>%
    group_by(MONTH) %>%
    summarise(
        travel_infectious = sum(travel_infectious),
        travel_infectious_lower = sum(travel_infectious_lower),
        travel_infectious_upper = sum(travel_infectious_upper),
        population = mean(population)
    ) %>%
    ungroup()

air_travel_from_nola_risk <- air_travel_from_nola %>%
    left_join(nola_cases, by=c("Date" = "MONTH")) %>%
    ungroup() %>%
    mutate(
        infectious_travellers = (travel_infectious/population),
        import_risk = infectious_travellers * n,
        infectious_travellers_lower = (travel_infectious_lower/population),
        import_risk = infectious_travellers_lower * n,
        infectious_travellers_upper = (travel_infectious_upper/population),
        import_risk = infectious_travellers_upper * n,
        Date = fct_relevel(as.factor(Date), c("February-2020", "March-2020", "April-2020", "May-2020"))
    ) %>%
    drop_na(`Destination Province`, import_risk)

air_travel_from_nola_risk %>%
    ggplot(aes(Date, reorder(`Destination Province`, -import_risk), fill = import_risk)) + geom_tile() + geom_text(aes(Date, reorder(`Destination Province`, -import_risk), label = round(import_risk, 2))) + scale_fill_distiller(palette = "Blues", direction = 1, trans="log10") + theme_bw()
ggsave("./plots/air_travel_risk_from_nola.svg", w= 15, h = 10)

air_travel_from_nola_risk %>%
    ggplot(aes(Date, reorder(group, -n), fill = n)) + geom_tile() + geom_text(aes(MONTH, reorder(group, -n), label = round(n, 2))) + scale_fill_distiller(palette = "Blues", direction = 1, trans="log10") + theme_bw()
ggsave("./plots/air_travel_from_nola.svg", w= 15, h = 10)
