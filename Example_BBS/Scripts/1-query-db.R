####################
# 1 - query database for bbs data
#
# REQUIRES PASSWORD - USE PROVIDED .rds OBJECT TO RECRATE DEPENDENT ANALYSES
# 
# Creates:
# * bbs_query_5930_2018-10-30.rds -> bbs results from query
####################


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)

# set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Example_BBS/Data/'))


# access DB ---------------------------------------------------------------

pass <- readLines('db_pass.txt')

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg, 
                      user = "cyoungflesh", 
                      password = pass, 
                      host = "35.221.16.125", 
                      port = 5432, 
                      dbname = "sightings")



# Get all routes run in 2016 ----------------------------------------------------

# All unique surveys
# Just 2016
# Only rpid = 101 and routetypedetailid = 1
# see metadata here: https://github.com/phenomismatch/sightings-database/blob/master/docs/metadata/bbs/RunType.pdf
# need to zero fill as 0 counts are meaningful


data <- DBI::dbGetQuery(cxn, paste0("SELECT DISTINCT ON (event_id) event_id, year, day, 
                                    lng, lat,
                                    (place_json ->> 'countrynum')::int AS country_number,
                                    (place_json ->> 'statenum')::int AS state_number,
                                    (place_json ->> 'route')::int AS route_number,
                                    (place_json ->> 'routename') AS route_name,
                                    (event_json ->> 'rpid')::int AS rpid,
                                    (event_json ->> 'runtype')::int AS runtype,
                                    (place_json ->> 'routetypedetailid')::int AS routetypedetailid
                                    FROM places
                                    JOIN events USING (place_id)
                                    JOIN counts USING (event_id)
                                    JOIN taxons USING (taxon_id)
                                    WHERE dataset_id = 'bbs'
                                    AND (place_json ->> 'routetypedetailid')::int = 1
                                    AND (event_json ->> 'rpid')::int = 101
                                    AND year = 2016;
                                    "))


#add species columns
data[c('target_count')] <- NA


#Northern Cardinal AOU
AOU <- '5930'


#query just cardinal
temp <- DBI::dbGetQuery(cxn, paste0("SELECT event_id, year, day, 
                                    lng, lat, count, common_name,
                                    (place_json ->> 'countrynum')::int AS country_number,
                                    (place_json ->> 'statenum')::int AS state_number,
                                    (place_json ->> 'route')::int AS route_number,
                                    (place_json ->> 'routename') AS route_name,
                                    (event_json ->> 'rpid')::int AS rpid,
                                    (event_json ->> 'runtype')::int AS runtype,
                                    (count_json ->> 'aou') AS aou,
                                    (place_json ->> 'routetypedetailid')::int AS routetypedetailid
                                    FROM places
                                    JOIN events USING (place_id)
                                    JOIN counts USING (event_id)
                                    JOIN taxons USING (taxon_id)
                                    WHERE dataset_id = 'bbs'
                                    AND (place_json ->> 'routetypedetailid')::int = 1
                                    AND (event_json ->> 'rpid')::int = 101
                                    AND (count_json ->> 'aou') ~ '\\y",AOU,"\\y'
                                    AND year = 2016;
                                    "))

#indices in data that match species ind
ind <- which(data$event_id %in% temp$event_id)

data[ind,'target_count'] <- temp$count

#indices to fill with 0s (observations of species not made for these events)
n_ind <- (1:NROW(data))[-ind]

data[n_ind,'target_count'] <- 0

query_rds_name <- paste0('bbs_query_', AOU, '_', Sys.Date(), '.rds')

saveRDS(data, paste0(query_rds_name))
