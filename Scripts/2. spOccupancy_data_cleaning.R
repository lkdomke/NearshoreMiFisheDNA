# Lia Domke
# 5-14-25
# Clean data so that it nicely can be formatted and input into species occupancy models 


#long <- read.csv("Data/combined_gear_wcomm_long_4-2-25.csv")
# this has the site info for the paired sites only
site.info2 <- read.csv("Data/combined_metadata_4-2-25.csv")

# pull in the long form data with all sites (35 paired + 2 edna only)
long.full <- read.csv("Data/combined_gear_wcomm_extra-eDNA_5-9-25.csv") %>%
  dplyr::select(-X)

# get the lat/lon for the edna only sites in 2021
site <- read.csv("Data/Site_metadata_allyears_8-2022.csv") %>%
  unite(bay_id, c(bay_code:bay_sample)) %>%
  filter(bay_id %in% c("BTRB_A", "CLAM_A", "KLWA_A", "NFEI_B")) %>%
  dplyr::select(place_name, bay_id, latitude, longitude) %>%
  mutate(date = ifelse(bay_id == "BTRB_A", "2021-08-07", NA),
         date = ifelse(bay_id == "CLAM_A", "2021-08-08", date),
         date = ifelse(bay_id == "KLWA_A", "2021-08-06", date),
         date = ifelse(bay_id == "NFEI_B", "2021-08-09", date)) %>%
  mutate(date = ymd(date),
         year = year(date),
         julian = yday(date)) %>%
  unite(bay_year2, c(bay_id, year), remove = F)

# site has the additional lat/lon for the edna sites, but doesn't have the doy 
edna.only <- filter(long.full, bay_year2 %in% c("BTRB_A_2021","NFEI_B_2021" #, "CLAM_A_2021", "KLWA_A_2021"
)) %>%
  distinct(bay_year, bay_year2, intx, habitat, type)

# add in the site info
edna.only.sites <- edna.only %>%
  left_join(site) %>%
  dplyr::select(-bay_id)

edna.only.fish <- long.full %>%
  filter(bay_year2 %in% c("BTRB_A_2021", "NFEI_B_2021" #, "CLAM_A_2021", "KLWA_A_2021" 
  )) %>%
  mutate(date = ymd(date)) %>%
  dplyr::select(-c(place_name, latitude, longitude, date, year, julian)) %>%
  full_join(edna.only.sites, by = c("bay_year", "bay_year2", "intx", "habitat", "type"))

long.full2 <- long.full %>%
  # removes all edna sites - well add back in btrb and nfei but leave out clam/klwa 2021
  filter(!(bay_year2 %in% c("BTRB_A_2021","NFEI_B_2021", "CLAM_A_2021", "KLWA_A_2021" 
  ))) %>%
  mutate(date = ymd(date))

# go back and create the clean long dataframe with species info and site info 
long <- rbind(long.full2, edna.only.fish) %>%
  dplyr::select(-c(place_name, date))

distinct(long, bay_year2, habitat, julian) # 37 sites with edna or seine (35 paired)

# create a site metadata
sites <- long %>%
  distinct(bay_year2, habitat, julian, latitude, longitude,year)

## write out clean data -- 
# write.csv(long, "Data/spOccupancy_fishdat_long_5-14-25.csv")
# write.csv(sites, "Data/spOccupancy_fishdat_sites_5-14-25.csv")
