
#Social relationships promote access to food and information in wild jackdaws  
#Proceedings of the Royal Society B: Biological Sciences
#10.1098/rspb. 2026-1462

#Author: Luca Hahn
#Last update: 26/06/2026

#(1) DATA PREPARATION ----

#Load packages ----

library(brms)
library(car)
library(carData)
library(chisq.posthoc.test)
library(ClusterR)
library(corrplot)
library(data.table)
library(DHARMa)
library(dplyr)
library(emmeans)
library(extrafont)
font_import()
loadfonts()
library(ggeffects)
library(ggplot2)
library(ggmap)
library(ggpubr)
library(ggrepel)
library(ggsn)
library(ggthemes)
library(gridGraphics)
library(grid)
library(gridBase)
library(glmmTMB)
library(hms)
library(igraph)
library(interactions)
library(lme4)
library(lmerTest)
library(lubridate)
library(MASS)
library(modelr)
library(multcomp)
library(NBDA)
library(patchwork)
library(performance)
library(RColorBrewer)
library(reshape2)
library(rptR)
library(sjlabelled)
library(sjPlot)
library(stringi)
library(STbayes)
library(svMisc)
library(tidybayes)
library(tidyverse)

#Colourblind 
palette.colors(palette = "Okabe-Ito")

##         black        orange       skyblue   bluishgreen        yellow 
##     "#000000"     "#E69F00"     "#56B4E9"     "#009E73"     "#F0E442" 
##          blue    vermillion reddishpurple          gray 
##     "#0072B2"     "#D55E00"     "#CC79A7"     "#999999" 

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

#Read 2023 visit data ----

#Shows milliseconds
op <- options(digits.secs=1)

#Use directory where you want to look for RT files, concatenate paths
RT_files <- list.files("Data/2023", pattern = "RT", recursive = TRUE)
RT_paths <- paste("Data/2023",RT_files, sep = "/")

library(dplyr)
library(stringi)

# Unique days/positions
day_positions <- unique(substr(RT_files, 6, 16))
collapsed_list <- list()
visit_time <- 10  # seconds, gap to start new visit

for (i in seq_along(day_positions)) {
  
  progress(i, max.value = length(day_positions))
  
  day_files_idx <- which(substr(RT_files, 6, 16) == day_positions[i])
  
  # Read all files for this day
  day_position_list <- lapply(day_files_idx, function(j) {
    temp_day_file <- read.delim(RT_paths[j], header = TRUE, stringsAsFactors = FALSE)[,1:13]
    temp_day_file$feeder <- substr(RT_files[j], 6, 9)
    temp_day_file
  })
  
  temp_RT <- do.call(rbind, day_position_list)
  
  # Ensure Hmsec is numeric, replace NA with 0
  temp_RT$Hmsec <- as.numeric(temp_RT$Hmsec)
  temp_RT$Hmsec[is.na(temp_RT$Hmsec)] <- 0
  
  # Pad 2-digit Hmsec by multiplying by 10
  pad_idx <- which(temp_RT$Hmsec >= 10 & temp_RT$Hmsec < 100)
  temp_RT$Hmsec[pad_idx] <- temp_RT$Hmsec[pad_idx] * 10
  
  # Build exact timestamps including milliseconds
  temp_RT$Time_exact <- as.POSIXct(temp_RT$Date, format="%Y-%m-%d %H:%M:%S", tz="UTC") +
    (temp_RT$Hmsec / 1024)
  
  # Keep only valid tags
  temp_tags <- temp_RT %>% filter(nchar(TagID_hex) == 10)
  if (nrow(temp_tags) == 0) next
  
  # Create read intervals with 0.5*(Reps - 1) adjustment
  reads <- data.frame(
    Time  = temp_tags$Time_exact,
    End   = temp_tags$Time_exact + 0.5 * (temp_tags$Reps - 1),
    tag   = temp_tags$TagID_hex,
    feeder= temp_tags$feeder,
    event = temp_tags$Event
  )
  
  reads <- reads %>% arrange(feeder, Time)
  
  visits_all <- list()
  visit_counter <- 1
  
  # Loop over feeders
  for (f in unique(reads$feeder)) {
    feeder_reads <- reads %>% filter(feeder == f) %>% arrange(Time)
    if (nrow(feeder_reads) == 0) next
    
    current_tag   <- feeder_reads$tag[1]
    current_start <- feeder_reads$Time[1]
    current_end   <- feeder_reads$End[1]
    current_event <- feeder_reads$event[1]
    
    for (r in 2:nrow(feeder_reads)) {
      this_tag   <- feeder_reads$tag[r]
      this_time  <- feeder_reads$Time[r]
      this_end   <- feeder_reads$End[r]
      this_event <- feeder_reads$event[r]
      
      if (any(is.na(c(this_tag, this_time, this_end)))) next
      
      gap <- as.numeric(difftime(this_time, current_end, units = "secs"))
      
      # New visit condition
      if (this_tag != current_tag || gap > visit_time) {
        
        # Save previous visit
        visits_all[[visit_counter]] <- data.frame(
          feeder = f,
          tag    = current_tag,
          start  = current_start,
          end    = current_end,
          event  = current_event
        )
        visit_counter <- visit_counter + 1
        
        # Start new visit
        current_tag   <- this_tag
        current_start <- this_time
        current_end   <- this_end
        current_event <- this_event
        
      } else {
        # Continue same visit
        current_end <- max(current_end, this_end)
      }
    }
    
    # Save last visit for feeder
    visits_all[[visit_counter]] <- data.frame(
      feeder = f,
      tag    = current_tag,
      start  = current_start,
      end    = current_end,
      event  = current_event
    )
    visit_counter <- visit_counter + 1
  }
  
  collapsed_list[[i]] <- dplyr::bind_rows(visits_all)
}

# Combine all visits
visit_data <- dplyr::bind_rows(collapsed_list)

setDT(visit_data)

# Ensure timestamps are POSIXct with milliseconds
visit_data[, start := ymd_hms(start)]
visit_data[, end   := ymd_hms(end)]
options(digits.secs = 3)

# Assign a unique event ID per visit
visit_data[, event_id := paste(feeder, format(start, "%Y-%m-%d %H:%M:%OS3"), format(end, "%Y-%m-%d %H:%M:%OS3"), tag, sep = " ")]

#Identify duplicates
dup_rows <- visit_data %>%
  group_by(event_id) %>%
  filter(n() > 1) %>%
  ungroup()

#Remove duplicates
visit_data <- visit_data %>%
  distinct(event_id, .keep_all = TRUE)

# Ensure milliseconds are shown
options(digits.secs = 3)

# Check result
head(visit_data$start)
head(visit_data$end)

#Fix remaining instances of visit overlap
visit_data <- visit_data %>%
  arrange(feeder, start) %>%
  group_by(feeder) %>%
  mutate(
    next_start = lead(start),
    end = dplyr::if_else(
      !is.na(next_start) & end > next_start,
      next_start,
      as.POSIXct(end, origin = "1970-01-01")
    )
  ) %>%
  dplyr::select(-next_start) %>%
  ungroup()

# Optional: compute time differences between consecutive visits per feeder
visit_data <- visit_data %>%
  arrange(feeder, start) %>%
  group_by(feeder) %>%
  mutate(time_diff_pre = as.numeric(start - lag(end))) %>%
  ungroup()

sum(visit_data$time_diff_pre < 0, na.rm = T)

#Find instances in which time is still NA
visit_data_NA <- subset(visit_data, is.na(visit_data$start))

#Filter instances in which time is not NA
visit_data <- subset(visit_data, !is.na(visit_data$start))

#Remove instances where JID = NA and test tags
visit_data$JID <- LH_RFID$ID[match(visit_data$tag, LH_RFID$RFID)]
visit_data <- subset(visit_data, JID != "NA")

visit_data23 <- visit_data 
visit_data23$year <- 2023

#Read 2024 visit data ----

#Shows milliseconds
op <- options(digits.secs=3)

#Use directory where you want to look for RT files, concatenate paths
RT_files <- list.files("Data/2024", pattern = "RT", recursive = TRUE)
RT_paths <- paste("Data/2024",RT_files, sep = "/")

library(dplyr)
library(stringi)

# Unique days/positions
day_positions <- unique(substr(RT_files, 6, 16))
collapsed_list <- list()
visit_time <- 10  # seconds, gap to start new visit

for (i in seq_along(day_positions)) {
  
  progress(i, max.value = length(day_positions))
  
  day_files_idx <- which(substr(RT_files, 6, 16) == day_positions[i])
  
  # Read all files for this day
  day_position_list <- lapply(day_files_idx, function(j) {
    temp_day_file <- read.delim(RT_paths[j], header = TRUE, stringsAsFactors = FALSE)[,1:13]
    temp_day_file$feeder <- substr(RT_files[j], 6, 9)
    temp_day_file
  })
  
  temp_RT <- do.call(rbind, day_position_list)
  
  # Ensure Hmsec is numeric, replace NA with 0
  temp_RT$Hmsec <- as.numeric(temp_RT$Hmsec)
  temp_RT$Hmsec[is.na(temp_RT$Hmsec)] <- 0
  
  # Pad 2-digit Hmsec by multiplying by 10
  pad_idx <- which(temp_RT$Hmsec >= 10 & temp_RT$Hmsec < 100)
  temp_RT$Hmsec[pad_idx] <- temp_RT$Hmsec[pad_idx] * 10
  
  # Build exact timestamps including milliseconds
  temp_RT$Time_exact <- as.POSIXct(temp_RT$Date, format="%Y-%m-%d %H:%M:%S", tz="UTC") +
    (temp_RT$Hmsec / 1024)
  
  # Keep only valid tags
  temp_tags <- temp_RT %>% filter(nchar(TagID_hex) == 10)
  if (nrow(temp_tags) == 0) next
  
  # Create read intervals with 0.5*(Reps - 1) adjustment
  reads <- data.frame(
    Time  = temp_tags$Time_exact,
    End   = temp_tags$Time_exact + 0.5 * (temp_tags$Reps - 1),
    tag   = temp_tags$TagID_hex,
    feeder= temp_tags$feeder,
    event = temp_tags$Event
  )
  
  reads <- reads %>% arrange(feeder, Time)
  
  visits_all <- list()
  visit_counter <- 1
  
  # Loop over feeders
  for (f in unique(reads$feeder)) {
    feeder_reads <- reads %>% filter(feeder == f) %>% arrange(Time)
    if (nrow(feeder_reads) == 0) next
    
    current_tag   <- feeder_reads$tag[1]
    current_start <- feeder_reads$Time[1]
    current_end   <- feeder_reads$End[1]
    current_event <- feeder_reads$event[1]
    
    for (r in 2:nrow(feeder_reads)) {
      this_tag   <- feeder_reads$tag[r]
      this_time  <- feeder_reads$Time[r]
      this_end   <- feeder_reads$End[r]
      this_event <- feeder_reads$event[r]
      
      if (any(is.na(c(this_tag, this_time, this_end)))) next
      
      gap <- as.numeric(difftime(this_time, current_end, units = "secs"))
      
      # New visit condition
      if (this_tag != current_tag || gap > visit_time) {
        
        # Save previous visit
        visits_all[[visit_counter]] <- data.frame(
          feeder = f,
          tag    = current_tag,
          start  = current_start,
          end    = current_end,
          event  = current_event
        )
        visit_counter <- visit_counter + 1
        
        # Start new visit
        current_tag   <- this_tag
        current_start <- this_time
        current_end   <- this_end
        current_event <- this_event
        
      } else {
        # Continue same visit
        current_end <- max(current_end, this_end)
      }
    }
    
    # Save last visit for feeder
    visits_all[[visit_counter]] <- data.frame(
      feeder = f,
      tag    = current_tag,
      start  = current_start,
      end    = current_end,
      event  = current_event
    )
    visit_counter <- visit_counter + 1
  }
  
  collapsed_list[[i]] <- dplyr::bind_rows(visits_all)
}

# Combine all visits
visit_data <- dplyr::bind_rows(collapsed_list)

setDT(visit_data)

# Ensure timestamps are POSIXct with milliseconds
visit_data[, start := ymd_hms(start)]
visit_data[, end   := ymd_hms(end)]
options(digits.secs = 3)

# Assign a unique event ID per visit
visit_data[, event_id := paste(feeder, format(start, "%Y-%m-%d %H:%M:%OS3"), format(end, "%Y-%m-%d %H:%M:%OS3"), tag, sep = " ")]

#Identify duplicates
dup_rows <- visit_data %>%
  group_by(event_id) %>%
  filter(n() > 1) %>%
  ungroup()

#Remove duplicates
visit_data <- visit_data %>%
  distinct(event_id, .keep_all = TRUE)

# Ensure milliseconds are shown
options(digits.secs = 3)

# Check result
head(visit_data$start)
head(visit_data$end)

#Fix remaining instances of visit overlap
visit_data <- visit_data %>%
  arrange(feeder, start) %>%
  group_by(feeder) %>%
  mutate(
    next_start = lead(start),
    end = dplyr::if_else(
      !is.na(next_start) & end > next_start,
      next_start,
      as.POSIXct(end, origin = "1970-01-01")
    )
  ) %>%
  dplyr::select(-next_start) %>%
  ungroup()

# Optional: compute time differences between consecutive visits per feeder
visit_data <- visit_data %>%
  arrange(feeder, start) %>%
  group_by(feeder) %>%
  mutate(time_diff_pre = as.numeric(start - lag(end))) %>%
  ungroup()

sum(visit_data$time_diff_pre < 0, na.rm = T)

#Find instances in which time is still NA
visit_data_NA <- subset(visit_data, is.na(visit_data$start))

#Filter instances in which time is not NA
visit_data <- subset(visit_data, !is.na(visit_data$start))

#Remove instances where JID = NA and test tags
visit_data$JID <- LH_RFID$ID[match(visit_data$tag, LH_RFID$RFID)]

visit_data <- subset(visit_data, JID != "NA")

visit_data24 <- visit_data 
visit_data24$year <- 2024

#Combine 2023 and 2024 data 
visit_data <- rbind(visit_data23, visit_data24) 

# Optional: compute time differences between consecutive visits per feeder
visit_data <- visit_data %>%
  arrange(feeder, start) %>%
  group_by(feeder) %>%
  mutate(time_diff_pre = as.numeric(start - lag(end))) %>%
  ungroup()

sum(visit_data$time_diff_pre < 0, na.rm = T)

visit_data$array <- NULL

#Adding information about visit duration
visit_data$interval <- interval(visit_data$start,visit_data$end)
visit_data$visit_duration <- as.duration(visit_data$interval)

#Adding information about perch: primary or secondary
visit_data[substr(visit_data$feeder,4,4)=="1","perch"]<-"primary"
visit_data[substr(visit_data$feeder,4,4)=="2","perch"]<-"secondary"

#Adding information about site: X, Y, Z
visit_data[substr(visit_data$feeder,1,1)=="Y","site"]<-"Y"
visit_data[substr(visit_data$feeder,1,1)=="Z","site"]<-"Z"

#Adding information about feeder position 
visit_data[substr(visit_data$feeder,1,2)=="Y1","position"]<-"Y1"
visit_data[substr(visit_data$feeder,1,2)=="Y2","position"]<-"Y2"
visit_data[substr(visit_data$feeder,1,2)=="Y3","position"]<-"Y3"
visit_data[substr(visit_data$feeder,1,2)=="Y4","position"]<-"Y4"
visit_data[substr(visit_data$feeder,1,2)=="Y5","position"]<-"Y5"
visit_data[substr(visit_data$feeder,1,2)=="Y6","position"]<-"Y6"
visit_data[substr(visit_data$feeder,1,2)=="Z1","position"]<-"Z1"
visit_data[substr(visit_data$feeder,1,2)=="Z2","position"]<-"Z2"
visit_data[substr(visit_data$feeder,1,2)=="Z3","position"]<-"Z3"
visit_data[substr(visit_data$feeder,1,2)=="Z4","position"]<-"Z4"
visit_data[substr(visit_data$feeder,1,2)=="Z5","position"]<-"Z5"
visit_data[substr(visit_data$feeder,1,2)=="Z6","position"]<-"Z6"

#Adding day of year and day of study period 
visit_data$day <- yday(visit_data$start) #day of year
#visit_data$study_day <- visit_data$day - 77 #day of study period

#Remove single instance of time slip where day > 250
visit_data <- subset(visit_data, day < 250)

#Adding information about time of the day 
visit_data$time <- as_hms(visit_data$start)
visit_data$hour <- hour(visit_data$start)

#Remove single instance of time slip where time = 0
visit_data <- subset(visit_data, hour > 0)

#Remove time slip where year = 2033 
visit_data <- subset(visit_data, !substr(visit_data$start, 3, 4) == "33")

visit_data$partnerID <- LH_partner$PARTNER.ID[match(visit_data$JID,LH_partner$ID)]
visit_data$pairID <- LH_pairs$pair_ID[match(visit_data$JID, LH_pairs$ID)]
visit_data$motherID <- LH_mother$MOTHER.ID[match(visit_data$JID,LH_mother$ID)]
visit_data$fatherID <- LH_father$FATHER.ID[match(visit_data$JID,LH_father$ID)]

#Adding individual age
Current_year <- 2024  #Set reference year

LH$DATE = as.Date(LH$DATE, format="%d-%m-%Y")
LH$year <- as.numeric(format(LH$DATE,"%Y"))   #  #Get year from the date
Ringed <- LH[!LH$year == 2025,] %>% filter(CODE == "RINGED")   #Get only records of when birds were ringed for the first time

Ringed$known_age <- 0   #binary 0/1 do we know the exact age (e.g. birds ringed as a 6 are 0)
Ringed$min_age <- 0  #Either actual age (if known_age = 1), or minimum age (if known_age = 0) - currently as number of new years crossed.

for (i in  1:nrow(Ringed)) {
  if(Ringed[i,12] == "4"){
    Ringed$known_age[i] = 0
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +2)
  }
  else if(Ringed[i,12] %in% c("1","1J","3","3J")){
    Ringed$known_age[i] = 1
    Ringed$min_age[i] = (Current_year - Ringed$year[i]+1) 
  }
  else if(Ringed[i,12] == "5"){
    Ringed$known_age[i] = 1
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +2) 
  }
  else if (Ringed[i,12] == "6"){
    Ringed$known_age[i] = 0
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +3) 
  }
  else { Ringed$known_age[i] = 0
  Ringed$min_age[i] = NA }
}

sum(Ringed$known_age)
length(Ringed$known_age)

Ringed$age24 <- Ringed$min_age
Ringed$age23 <- Ringed$min_age - 1

unique_ID = unique(visits$JID)

visit_data$age <- ifelse(visit_data$year == 2024, Ringed$age24[match(visit_data$JID,Ringed$ID)], Ringed$age23[match(visit_data$JID,Ringed$ID)])
visit_data23$age <- Ringed$min_age[match(visit_data23$JID,Ringed$ID)] - 1
visit_data24$age <- Ringed$min_age[match(visit_data24$JID,Ringed$ID)]

#Read other datasets ----

#Load saved life history csv file 
LH <- read.csv("Data/LH20251207.csv", header = T, stringsAsFactors = F)

LH$DATE <- strptime(LH$DATE,format="%d/%m/%Y")
LH$DATE <- as.Date(LH$DATE, format = "%d/%m/%Y") # convert to date
LH$year <- as.numeric(format(LH$DATE,"%Y"))   

#Separate LH file with (1) pairs and then (2) adding pair ID
#(i) LH file with pairs and boxes
LH_pairs <- subset(LH[c("DATE", "ID", "SEX", "BOX", "PARTNER.ID", "year")])
LH_pairs <- LH_pairs[LH_pairs$PARTNER.ID != "",]
length(table(LH_pairs$ID)) #570 individuals with pair data

LH_pairs23 <- subset(LH_pairs, year < 2024) #all until 2023
LH_pairs24 <- subset(LH_pairs, year < 2025)
LH_pairs25 <- LH_pairs

#Get most recent entries for pair ID  (2023 and 2024)
LH_pairs23 <- LH_pairs23 %>%
  group_by(ID) %>% 
  arrange(desc(DATE)) %>% 
  slice(1:1)

LH_pairs24 <- LH_pairs24 %>%
  group_by(ID) %>% 
  arrange(desc(DATE)) %>% 
  slice(1:1)

LH_pairs25 <- LH_pairs25 %>%
  group_by(ID) %>% 
  arrange(desc(DATE)) %>% 
  slice(1:1)

#(ii) Adding pair ID
LH_pairs23  <- LH_pairs23 %>% 
  as_tibble() %>% 
  mutate(pair_ID = if_else(LH_pairs23$SEX == "F", paste(ID, PARTNER.ID), paste(PARTNER.ID, ID)))

LH_pairs24  <- LH_pairs24 %>% 
  as_tibble() %>% 
  mutate(pair_ID = if_else(LH_pairs24$SEX == "F", paste(ID, PARTNER.ID), paste(PARTNER.ID, ID)))

LH_pairs25  <- LH_pairs25 %>% 
  as_tibble() %>% 
  mutate(pair_ID = if_else(LH_pairs25$SEX == "F", paste(ID, PARTNER.ID), paste(PARTNER.ID, ID)))

#LH_box (info about box owners and juveniles born in boxes)
LH_box_owners <- LH[LH$BOX != "",]
LH_box_owners$year <- as.numeric(format(LH_box_owners$DATE,"%Y"))   
LH_box_owners$JID_year <- paste(LH_box_owners$ID, LH_box_owners$year, sep = "_")

LH_box_juv <- LH[LH$CODE == "FLEDGE",]
LH_box_juv$JID_year <- paste(LH_box_juv$ID, LH_box_juv$year, sep = "_")

LH_box <- rbind(LH_box_owners, LH_box_juv)

#LH_sex (summary with all individuals' sex)
LH_sex <- LH[LH$SEX != "",]

#LH body 
LH_body <- subset(LH, !LH$TARSUS == "")

#Adding individual age
Current_year <- 2024  #Set reference year

LH$DATE = as.Date(LH$DATE, format="%d-%m-%Y")
LH$year <- as.numeric(format(LH$DATE,"%Y"))   #  #Get year from the date
Ringed <- LH[!LH$year == 2025,] %>% filter(CODE == "RINGED")   #Get only records of when birds were ringed for the first time

Ringed$known_age <- 0   #binary 0/1 do we know the exact age (e.g. birds ringed as a 6 are 0)
Ringed$min_age <- 0  #Either actual age (if known_age = 1), or minimum age (if known_age = 0) - currently as number of new years crossed.

for (i in  1:nrow(Ringed)) {
  if(Ringed[i,12] == "4"){
    Ringed$known_age[i] = 0
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +2)
  }
  else if(Ringed[i,12] %in% c("1","1J","3","3J")){
    Ringed$known_age[i] = 1
    Ringed$min_age[i] = (Current_year - Ringed$year[i]+1) 
  }
  else if(Ringed[i,12] == "5"){
    Ringed$known_age[i] = 1
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +2) 
  }
  else if (Ringed[i,12] == "6"){
    Ringed$known_age[i] = 0
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +3) 
  }
  else { Ringed$known_age[i] = 0
  Ringed$min_age[i] = NA }
}

sum(Ringed$known_age)
length(Ringed$known_age)

Ringed$age24 <- Ringed$min_age
Ringed$age23 <- Ringed$min_age - 1

#LH_mother and LH_father
LH_mother <- LH[LH$MOTHER.ID != "",]
LH_father <- LH[LH$FATHER.ID != "",]

#LH_ring (to infer fledge dates for juveniles)
LH_ring <- LH[LH$CODE == "RINGED",]
LH_ring$day_ringed <- yday(LH_ring$DATE)

#LH_fledge (to infer fledge dates for juveniles)
LH_fledge <- LH[LH$CODE == "FLEDGE",]
LH_fledge$day_fledge <- yday(LH_fledge$DATE)

#Box distance (pairwise distances between boxes)
box_distance <- read.csv("Data/distance_matrix.csv", header = T, stringsAsFactors = F)
box_distance$box_dyad_ID <- ifelse(box_distance$InputID < box_distance$TargetID, paste(box_distance$InputID, box_distance$TargetID, sep = "_"), paste(box_distance$TargetID, box_distance$InputID, sep = "_"))

#Social pedigree (pairwise relatedness between birds)
social_pedigree <- read.csv("Data/social_pedigree.csv", header = T, stringsAsFactors = F)

#How many birds trapped in 2023 were already marked?
LH_captured <- subset(LH, LH$CODE == "RETRAP" | LH$CODE == "RINGED")
LH_captured$year <- substr(LH_captured$DATE, 1, 4)
LH_captured$JID_year <- paste(LH_captured$ID, LH_captured$year, sep = "_")
table(LH_captured$CODE[LH_captured$year ==2023], LH_captured$AGE[LH_captured$year ==2023])

#Read data on foraging associations obtained from Gaussian Mixture Model (to infer foraging associates from the year before)
edgesf22 <- read.csv("Data/edgesf22.csv", header = T, stringsAsFactors = F)
edgesf23 <- read.csv("Data/edgesf23.csv", header = T, stringsAsFactors = F)

edgesf22$dyad_ID <- ifelse(edgesf22$from < edgesf22$to, paste(edgesf22$from, edgesf22$to, sep = "_"), paste(edgesf22$to, edgesf22$from, sep = "_"))
edgesf23$dyad_ID <- ifelse(edgesf23$from < edgesf23$to, paste(edgesf23$from, edgesf23$to, sep = "_"), paste(edgesf23$to, edgesf23$from, sep = "_"))

#Visit data 
visit_data <- read.csv("Data/visit_data.csv", header = T, stringsAsFactors = F)

#Identify dyadic events ----

## 0. Load packages
library(data.table)
library(lubridate)

## 1. Prepare data

#Calculate visit interval
visit_data$interval <- interval(visit_data$start,visit_data$end)

#Take out all but "active" feeders (i.e. feeders with dual perches)
visit_data_dual <- arrange(visit_data %>% filter(feeder %in% c("Y1.1","Y1.2", "Y2.1", "Y2.2", "Y3.1", "Y3.2", "Y4.1", "Y4.2", "Z2.1", "Z2.2", "Z3.1", "Z3.2", "Z5.1", "Z5.2", "Z6.1", "Z6.2") & position %in% c("Y1", "Y2", "Y3", "Y4", "Z2", "Z3", "Z5", "Z6")),start)
visit_data_solo <- arrange(visit_data %>% filter(feeder %in% c("Y5.1", "Y6.1", "Z1.1", "Z4.1")),start)

#If all feeders are used
#visit_data

setDT(visit_data_dual)

#Ensure timestamps are POSIXct with milliseconds
visit_data_dual[, start := ymd_hms(start)]
visit_data_dual[, end   := ymd_hms(end)]
options(digits.secs = 3)

#Assign a unique event ID per visit
visit_data_dual[, event_id := paste(feeder, format(start, "%Y-%m-%d %H:%M:%OS3"), JID, sep = " ")]

#Remove invalid intervals
visit_data_dual <- visit_data_dual[end >= start]

#Sort for reproducibility
setorder(visit_data_dual, position, start, end)

# 2. Split into primary vs secondary perches
visit_data_primary   <- visit_data_dual[perch == "primary"]
visit_data_secondary <- visit_data_dual[perch == "secondary"]  # or feeder == "secondary"

# 3. Set keys for foverlaps
setkey(visit_data_primary, position, start, end)
setkey(visit_data_secondary, position, start, end)

# 4. Overlap join: secondary onto primary
hits <- foverlaps(visit_data_secondary, visit_data_primary, type="any", nomatch=0L)

# 5. Compute initiator / joiner and overlap
dual_events <- hits[tag != i.tag]

dual_events[, `:=`(
    initiator = fifelse(start <= i.start, tag, i.tag),
    joiner    = fifelse(start <= i.start, i.tag, tag),
    
    initiator_bout_start = fifelse(start <= i.start, start, i.start),
    initiator_bout_end   = fifelse(start <= i.start, end, i.end),
    
    joiner_bout_start = fifelse(start <= i.start, i.start, start),
    joiner_bout_end   = fifelse(start <= i.start, i.end, end),
    
    overlap_start = pmax(start, i.start),
    overlap_end   = pmin(end, i.end)
   )]

# Compute duration in seconds
dual_events[, duration := as.numeric(overlap_end - overlap_start, units = "secs")]

#Add a small constant to all overlaps (optional)
#constant_add <- 0.5  # seconds
#dual_events[, duration := duration + constant_add]
#dual_events[, overlap_end := overlap_start + duration]

dual_events <- dual_events[duration >= 0]

# 6. Assign primary perch event ID
# Primary perch is always the second table in the join (i)
dual_events[, primary_event_id := event_id]
dual_events[, secondary_event_id := i.event_id]

# 7. Build final dual_events dataset
dual_events <- dual_events[, .(
  initiator,
  joiner,
  position,
  initiator_feeder = ifelse(initiator == tag, feeder, i.feeder),
  joiner_feeder    = ifelse(joiner == tag, feeder, i.feeder),
  start = overlap_start,
  end   = overlap_end,
  duration,
  initiator_bout_start,
  initiator_bout_end,
  joiner_bout_start,
  joiner_bout_end,
  year,
  event_id = primary_event_id,
  secondary_event_id = secondary_event_id
)]


#Additional variables ----

#Add initiator and joiner perch
dual_events$initiator_perch <- ifelse(substr(dual_events$initiator_feeder, 4,4) == 1, "primary", "secondary")
dual_events$joiner_perch <- ifelse(substr(dual_events$joiner_feeder, 4,4) == 1, "primary", "secondary")

#Add initiator and joiner JID
dual_events$initiator_JID <- LH$ID[match(dual_events$initiator,LH$RFID)]
dual_events$joiner_JID <- LH$ID[match(dual_events$joiner,LH$RFID)]

#Add primary and secondary perch JID (jackdaw ID)
dual_events$primary_perch_JID <- ifelse(dual_events$initiator_perch == "primary", dual_events$initiator_JID, dual_events$joiner_JID)
dual_events$secondary_perch_JID <- ifelse(dual_events$initiator_perch == "primary", dual_events$joiner_JID, dual_events$initiator_JID)

#Add primary and secondary perch id (perch ID)
dual_events$primary_perch_id <- ifelse(dual_events$initiator_perch == "primary", dual_events$initiator_feeder, dual_events$joiner_feeder)
dual_events$secondary_perch_id <- ifelse(dual_events$initiator_perch == "primary", dual_events$joiner_feeder, dual_events$initiator_feeder)

#Add initiator and joiner duration
dual_events$initiator_duration <- as.numeric(dual_events$initiator_bout_end - dual_events$initiator_bout_start)
dual_events$joiner_duration <- as.numeric(dual_events$joiner_bout_end - dual_events$joiner_bout_start)

#Add difference between initiator and joiner arrival
dual_events$initiator_arriv <- dual_events$initiator_bout_start
dual_events$joiner_arriv <- dual_events$joiner_bout_start
dual_events$arriv_diff <- as.numeric(dual_events$joiner_bout_start - dual_events$initiator_bout_start)

#Add difference between initiator and joiner departure
dual_events$initiator_dep <- dual_events$initiator_bout_end
dual_events$joiner_dep<- dual_events$joiner_bout_end
dual_events$dep_diff <- as.numeric(dual_events$joiner_dep - dual_events$initiator_dep)

#Add primary and secondary perch start
dual_events$primary_arriv <- ifelse(dual_events$initiator_perch == "primary", paste(dual_events$initiator_arriv), paste(dual_events$joiner_arriv))
dual_events$secondary_arriv <- ifelse(dual_events$initiator_perch == "primary", paste(dual_events$joiner_arriv), paste(dual_events$initiator_arriv))
dual_events$primary_duration <- ifelse(dual_events$initiator_perch == "primary", dual_events$initiator_duration, dual_events$joiner_duration)

#Add initiator-joiner ID (considers both variants per dyad ID)
dual_events$initiator_joiner_ID <- paste(dual_events$initiator_JID, dual_events$joiner_JID)

#Add unique dyad ID
dual_events$dyad_ID <- ifelse(dual_events$initiator_JID < dual_events$joiner_JID, 
                              paste(dual_events$initiator_JID, dual_events$joiner_JID, sep = "_"),
                              paste(dual_events$joiner_JID, dual_events$initiator_JID, sep = "_"))

#Add initiator and joiner partner ID
dual_events$initiator_partnerID23 <- LH_pairs23$PARTNER.ID[match(dual_events$initiator_JID,LH_pairs23$ID)]
dual_events$joiner_partnerID23 <- LH_pairs23$PARTNER.ID[match(dual_events$joiner_JID,LH_pairs23$ID)]

dual_events$initiator_partnerID24 <- LH_pairs24$PARTNER.ID[match(dual_events$initiator_JID,LH_pairs24$ID)]
dual_events$joiner_partnerID24 <- LH_pairs24$PARTNER.ID[match(dual_events$joiner_JID,LH_pairs24$ID)]

dual_events$initiator_partnerID25 <- LH_pairs25$PARTNER.ID[match(dual_events$initiator_JID,LH_pairs25$ID)]
dual_events$joiner_partnerID25 <- LH_pairs25$PARTNER.ID[match(dual_events$joiner_JID,LH_pairs25$ID)]

#Add pair status
dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(pair23 = if_else(joiner_JID == initiator_partnerID23 & initiator_JID == joiner_partnerID23,"pair", "non-pair"))

dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(pair24 = if_else(joiner_JID == initiator_partnerID24 & initiator_JID == joiner_partnerID24,"pair", "non-pair"))

dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(pair25 = if_else(joiner_JID == initiator_partnerID25 & initiator_JID == joiner_partnerID25,"pair", "non-pair"))

dual_events$pair <- ifelse(dual_events$year == 2023, dual_events$pair23, dual_events$pair24)

#Manually add some pairs that were visited together many times and were found to be paired in 2025
dual_events$pair[dual_events$dyad_ID == "J4587_J4634"] <- "pair"
dual_events$pair[dual_events$dyad_ID == "J4646_J5053"] <- "pair"

#Add initiator/joiner and year
dual_events$initiator_year <- paste(dual_events$initiator_JID, dual_events$year, sep = "_") 
dual_events$joiner_year <- paste(dual_events$joiner_JID, dual_events$year, sep = "_") 

#Add initiator/joiner sex
dual_events$initiator_sex <- LH_sex$SEX[match(dual_events$initiator,LH_sex$RFID)]
dual_events$joiner_sex <- LH_sex$SEX[match(dual_events$joiner,LH_sex$RFID)]

#Add sex combinations initiator and joiner (both complete and 'simple')
dual_events <- dual_events %>%
  mutate(sex_combination=case_when(
    initiator_sex =="M" & joiner_sex=="M" ~ "MM",
    initiator_sex =="M" & joiner_sex== "F" ~ "MF",
    initiator_sex =="F" & joiner_sex== "M" ~ "FM",
    initiator_sex =="F" & joiner_sex== "F" ~ "FF",
    initiator_sex =="F" & is.na(joiner_sex) ~"FNA",
    is.na(initiator_sex) & joiner_sex== "F" ~ "NAF",
    initiator_sex =="M" & is.na(joiner_sex) ~"MNA",
    is.na(initiator_sex) & joiner_sex== "M" ~ "NAM",
    is.na(initiator_sex) & is.na(joiner_sex) ~ "NANA"))

#Sex class combination simple: regardless of arrival order, just 
dual_events <- dual_events %>%
  mutate(sex_combination_simple =case_when(
    sex_combination =="FM" | sex_combination =="MF" ~ "FM",
    sex_combination =="FF" | sex_combination== "FF" ~ "FF",
    sex_combination =="MM" | sex_combination== "MM" ~ "MM",
    sex_combination =="MNA" | sex_combination== "NAM" ~ "MNA",
    sex_combination =="FNA" | sex_combination== "NAF" ~ "FNA",
    sex_combination =="NANA" | sex_combination== "NANA" ~ "NANA"))

#Add initiator and joiner age, and age difference  
dual_events$initiator_age <- Ringed$min_age[match(dual_events$initiator_JID, Ringed$ID)]
dual_events$joiner_age <- Ringed$min_age[match(dual_events$joiner_JID, Ringed$ID)]

dual_events$initiator_age <- ifelse(dual_events$year == 2024, dual_events$initiator_age, dual_events$initiator_age - 1)
dual_events$joiner_age <- ifelse(dual_events$year == 2024, dual_events$joiner_age, dual_events$joiner_age - 1)

dual_events$age_diff <- dual_events$initiator_age - dual_events$joiner_age
dual_events$age_diff_abs <- abs(dual_events$age_diff)

#Add age class: adult and juvenile
dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(initiator_age_class = if_else(initiator_age > 1,"adult", "juvenile"))

dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(joiner_age_class = if_else(joiner_age > 1,"adult", "juvenile"))

#Add age class combinations initiator and joiner
dual_events <- dual_events %>%
  mutate(age_class_combination=case_when(
    initiator_age_class =="adult" & joiner_age_class =="adult" ~ "AA",
    initiator_age_class =="adult" & joiner_age_class =="juvenile" ~ "AJ",
    initiator_age_class =="juvenile" & joiner_age_class =="adult" ~ "JA",
    initiator_age_class =="juvenile" & joiner_age_class =="juvenile" ~ "JJ"))

#Add age class combination simple (regardless of arrival time)
dual_events <- dual_events %>%
mutate(age_class_combination_simple =case_when(
  age_class_combination =="AA" ~ "AA",
  age_class_combination =="JA" ~ "AJ",
  age_class_combination =="AJ" ~ "AJ",
  age_class_combination =="JJ" ~ "JJ"))
  
#Add initiator and joiner mother ID
dual_events$initiator_motherID <- LH_mother$MOTHER.ID[match(dual_events$initiator,LH_mother$RFID)]
dual_events$joiner_motherID <- LH_mother$MOTHER.ID[match(dual_events$joiner,LH_mother$RFID)]

#Add initiator and joiner father ID
dual_events$initiator_fatherID <- LH_father$FATHER.ID[match(dual_events$initiator,LH_father$RFID)]
dual_events$joiner_fatherID <- LH_father$FATHER.ID[match(dual_events$joiner,LH_father$RFID)]

#Add offspring-parent
dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(offspring_parent = if_else(joiner_JID == initiator_motherID | joiner_JID == initiator_fatherID,"offspring_parent", "no"))

#Add parent-offspring
dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(parent_offspring = if_else(initiator_JID == joiner_motherID | initiator_JID == joiner_fatherID,"parent_offspring", "no"))

#Add kin (parent and offspring only, no siblings)
dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(parent_offspring_kin = if_else(parent_offspring == "parent_offspring" | offspring_parent == "offspring_parent","kin", "non-kin"))

#Add mother sibling
dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(mother_sibling = if_else(is.na(initiator_motherID) & is.na(joiner_motherID) & initiator_motherID == joiner_motherID , "sibling", "non-sibling"))

dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(mother_sibling = if_else(initiator_motherID == joiner_motherID , "sibling", "non-sibling"))

#Add father sibling
dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(father_sibling = if_else(is.na(initiator_fatherID) & is.na(joiner_fatherID) & initiator_fatherID == joiner_fatherID , "sibling", "non-sibling"))

dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(father_sibling = if_else(initiator_fatherID == joiner_fatherID , "sibling", "non-sibling"))

#Add sibling
dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(sibling = if_else(is.na(initiator_fatherID) & is.na(joiner_fatherID) & initiator_fatherID == joiner_fatherID , "sibling", "non-sibling"))

dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(sibling = if_else(mother_sibling == "sibling" | father_sibling == "sibling", "sibling", "non-sibling"))

#Add first-order kin overall 
dual_events <- dual_events %>% 
  as_tibble() %>% 
  mutate(kin = if_else(sibling == "sibling" | parent_offspring == "parent_offspring" | offspring_parent == "offspring_parent","kin", "non-kin"))

#Add nestbox ids
dual_events$initiator_box <- LH_box$BOX[match(dual_events$initiator_year,LH_box$JID_year)]
dual_events$joiner_box <- LH_box$BOX[match(dual_events$joiner_year,LH_box$JID_year)]
dual_events$joiner_box[dual_events$joiner_box == "X32"] <- "Y32" #correct an error

#Add binary box ownership
dual_events$initiator_box_binary <- ifelse(!is.na(dual_events$initiator_box), 1, 0)
dual_events$joiner_box_binary <- ifelse(!is.na(dual_events$joiner_box), 1, 0)

#Add box distance
dual_events$box_dyad_ID <- ifelse(dual_events$initiator_box < dual_events$joiner_box, paste(dual_events$initiator_box, dual_events$joiner_box, sep = "_"), paste(dual_events$joiner_box, dual_events$initiator_box, sep = "_"))
dual_events$box_distance <- box_distance$Distance[match(dual_events$box_dyad_ID, box_distance$box_dyad_ID)]

#Add neighbour 
dual_events$neighbour <- ifelse(dual_events$box_distance < 50, 1, 0)
dual_events$neighbour[is.na(dual_events$neighbour)] <- 0
table(dual_events$neighbour)

#Add site 
dual_events$site <- substr(dual_events$position, 1, 1)

#Add individuals' preferred foraging site
#Visits per individual per site and per feeder, preferred sites
perindivpersite <- visit_data  %>%  count(JID, site)

perindivprefsite <- perindivpersite %>% 
  group_by(JID) %>%
  slice(which.max(n))

dual_events$initiator_pref_site <- perindivprefsite$site[match(dual_events$initiator_JID, perindivprefsite$JID)]
dual_events$joiner_pref_site <- perindivprefsite$site[match(dual_events$joiner_JID, perindivprefsite$JID)]
dual_events$site_comb <- paste(dual_events$initiator_pref_site, dual_events$joiner_pref_site, sep = "")
dual_events$site_comb_simple <- ifelse(dual_events$site_comb == "YZ" | dual_events$site_comb == "ZY", "YZ", dual_events$site_comb)
dual_events$site_comb_binary <- ifelse(dual_events$initiator_pref_site == dual_events$joiner_pref_site, 1, 0)
table(dual_events$site_comb_binary)

#Add relationship type
dual_events <- dual_events %>%
  mutate(relationship =case_when(
    pair =="pair" & kin == "non-kin" ~ "pair",
    pair == "non-pair" & kin == "kin" ~ "kin",
    pair == "pair" & is.na(kin) ~ "pair",
    is.na(pair) & kin == "kin" ~ "kin",
    pair == "non-pair" & kin == "non-kin" ~ "other",
    is.na(pair) & is.na(kin) ~ "other",
    is.na(pair) & kin == "non-kin"~ "other",
    pair == "non-pair" & is.na(kin) ~ "other",
    pair == "pair" & kin == "kin" ~ "pair"))

#Add social pedigree 
social_pedigree$dyad_ID <- ifelse(social_pedigree$Id1 < social_pedigree$Id2, paste(social_pedigree$Id1, social_pedigree$Id2, sep = "_"), paste(social_pedigree$Id2, social_pedigree$Id1, sep = "_"))
dual_events$social_pedigree <- social_pedigree$Id2_relatedness[match(dual_events$dyad_ID, social_pedigree$dyad_ID)]

table(dual_events$relationship, dual_events$social_pedigree)

#Add foraging association from previous year 
dual_events$prev_edgef <- ifelse(dual_events$year == 2023, edgesf22$weight[match(dual_events$dyad_ID, edgesf22$dyad_ID)], edgesf23$weight[match(dual_events$dyad_ID, edgesf23$dyad_ID)])
dual_events$prev_edgef_binary <- ifelse(is.na(dual_events$prev_edgef), 0, 1)

#Modify relationship type to include social pedigree and neighbours 
dual_events <- dual_events %>%
  mutate(relationship =case_when(
    relationship =="pair" ~ "pair",
    relationship == "kin" ~ "kin",
    neighbour == 1  ~ "neighbour",
    relationship == "other" & !is.na(social_pedigree) ~ "kin",
    relationship == "other" ~ "other",
    is.na(relationship) ~ "other"))

#Modify relationship type to include previous foraging associates and site residents
dual_events <- dual_events %>%
  mutate(relationship =case_when(
    relationship =="pair" ~ "pair",
    relationship == "kin" ~ "kin",
    relationship == "neighbour"  ~ "neighbour",
    relationship == "other" & !is.na(social_pedigree) ~ "kin",
    relationship == "other" & prev_edgef_binary == 1 & site_comb_binary == 1 ~ "foraging associate",
    relationship == "other" & prev_edgef_binary == 0 & site_comb_binary == 1 ~ "site resident",
    relationship == "other" & site_comb_binary == 0 ~ "other site"))

table(dual_events$relationship)

#Add perches
dual_events[substr(dual_events$initiator_feeder,4,4)=="1","initiator_perch"]<-"primary"
dual_events[substr(dual_events$initiator_feeder,4,4)=="2","initiator_perch"]<-"secondary"

#Add day 
dual_events$day <- yday(dual_events$start) #day of year
#dual_events$study_day <- dual_events$day - 77 #day of study period 

#Add proportion of overlap compared to primary feeder visit duration
dual_events$prop_dual_overlap <- (dual_events$duration -1)/ as.numeric(dual_events$initiator_duration)

#Summaries and subsets of data ----

#Dual events unique dyads (and when juveniles present, and only FM)
dual_events_dyads <- distinct(dual_events, dyad_ID, .keep_all = TRUE)

#Visit number dyads
dual_event_visits_dyads <- as.data.frame(table(dual_events$dyad_ID))
dual_event_visits_dyads$dyad_ID <- dual_event_visits_dyads$Var1
dual_event_visits_dyads$visit_number_dyad <- dual_event_visits_dyads$Freq
dual_event_visits_dyads <- subset(dual_event_visits_dyads, select = -c(Var1, Freq))

dual_events_dyads$visit_number_dyad <- dual_event_visits_dyads$visit_number_dyad[match(dual_events_dyads$dyad_ID, dual_event_visits_dyads$dyad_ID)]

#Subset initiator perch = feeder perch only
dual_events_primary <- subset(dual_events, initiator_perch == 'primary')

#Period in which juveniles present 
visit_data <- visit_data %>%
  dplyr::select(-interval)

visit_data %>%
  group_by(year) %>%
  filter(age == 1) %>%
  slice(1)

#juveniles first visit on 22-06 in each year 
#22-06 day 173 in 2023, and day 174 in 2024
dual_events_juv <- subset(dual_events, day > 173)
dual_events_nojuv <- subset(dual_events, day < 173)

#Dual events unique dyads (and when juveniles present, and only FM)
dual_events_dyads_juv <- distinct(dual_events_juv, dyad_ID, .keep_all = TRUE)

#Visit number dyads
dual_events_visit_dyads_juv <- as.data.frame(table(dual_events$dyad_ID))
dual_events_visit_dyads_juv$dyad_ID <- dual_events_visit_dyads_juv$Var1
dual_events_visit_dyads_juv$visit_number_dyad <- dual_events_visit_dyads_juv$Freq
dual_events_visit_dyads_juv <- subset(dual_events_visit_dyads_juv, select = -c(Var1, Freq))

dual_events_dyads_juv$visit_number_dyad <- dual_events_visit_dyads_juv$visit_number_dyad[match(dual_events_dyads_juv$dyad_ID, dual_events_visit_dyads_juv$dyad_ID)]

#Pairs only
dual_events_pair <- subset(dual_events, pair == 'pair')
dual_events_pair_juv <- subset(dual_events_juv, pair == 'pair')
dual_events_pair_primary <- subset(dual_events_primary, pair == 'pair')

length(table(dual_events_pair$dyad_ID)) #45 pairs
length(table(dual_events_pair_juv$dyad_ID))

#Kin only
dual_events_kin <- subset(dual_events, kin == 'kin')
dual_events_kin_juv <- subset(dual_events_juv, kin == 'kin')
dual_events_kin_primary <- subset(dual_events_primary, relationship == 'kin')

table(dual_events_kin$parent_offspring) #176 parent offspring, 
table(dual_events_kin$offspring_parent) #75 offspring parent, 
table(dual_events_kin$age_class_combination_simple) #57 AA, 225 AJ, 36 JJ 
length(table(dual_events_kin$dyad_ID)) #93 kin dyads

#Other only 
dual_events_other <- subset(dual_events, relationship == 'other')
dual_events_other_primary <- subset(dual_events_primary, relationship == 'other')
dual_events_other_juv <- subset(dual_events_juv, relationship == "other")

#Not others
dual_events_nonother <- subset(dual_events, !relationship == 'other')
dual_events_nonother_primary <- subset(dual_events_primary, !relationship == 'other')

#Attributes of all individuals participating in dual events
#For network attributes run code in (5) on networks first
dual_events_individuals <- dual_events[,c("start", "initiator_JID","joiner_JID")]
dual_events_individuals <- cbind(dual_events_individuals[1], stack(dual_events_individuals[2:3]))
dual_events_individuals <- dual_events_individuals[!duplicated(dual_events_individuals$values), ]
dual_events_individuals$JID <- dual_events_individuals$values
dual_events_individuals <- merge(x = dual_events_individuals, y = perindivpersite2, by = "JID", all = TRUE)
dual_events_individuals <- subset(dual_events_individuals, !is.na(dual_events_individuals$start))
dual_events_individuals <- subset(dual_events_individuals[,c(1,3)])
dual_events_individuals$age <- Ringed$min_age[match(dual_events_individuals$JID,Ringed$ID)]
dual_events_individuals$sex <- LH_sex$SEX[match(dual_events_individuals$JID,LH_sex$ID)]
dual_events_individuals$pair_ID <- LH_pairs$pair_ID[match(dual_events_individuals$JID,LH_pairs$ID)]
dual_events_individuals$pref_site <- perindivprefsite$site[match(dual_events_individuals$JID,perindivprefsite$JID)]

#Average duration of individual at primary perch during dual events
dual_event_duration <- as.data.frame(aggregate(dual_events_primary$initiator_duration, by = list(dual_events_primary$initiator_JID), FUN = "mean", na.rm = TRUE))
dual_event_duration$JID <- dual_event_duration$Group.1
dual_event_duration$visit_duration <- dual_event_duration$x
dual_event_duration <- subset(dual_event_duration, select = -c(Group.1, x))

#Number of primary feeder visits per individual 
dual_event_visits <- as.data.frame(table(dual_events_primary$initiator_JID))
dual_event_visits$JID <- dual_event_visits$Var1
dual_event_visits$visit_number <- dual_event_visits$Freq
dual_event_visits <- subset(dual_event_visits, select = -c(Var1, Freq))

#Full feeder visit sequence ----
 
#For data on juveniles following their parents

#Only dual feeding stations
visit_data_dual <- arrange(visit_data %>% filter(feeder %in% c("Y1.1","Y1.2", "Y2.1", "Y2.2", "Y3.1", "Y3.2", "Y4.1", "Y4.2", "Z2.1", "Z2.2", "Z3.1", "Z3.2", "Z5.1", "Z5.2", "Z6.1", "Z6.2") & position %in% c("Y1", "Y2", "Y3", "Y4", "Z2", "Z3", "Z5", "Z6")),start)

setDT(visit_data_dual)

# Ensure timestamps are POSIXct with milliseconds
visit_data_dual[, start := ymd_hms(start)]
visit_data_dual[, end   := ymd_hms(end)]
options(digits.secs = 3)

# Assign a unique event ID per visit
visit_data_dual[, event_id := paste(feeder, format(start, "%Y-%m-%d %H:%M:%OS3"), JID, sep = " ")]

#Visits at primary perches of dual feeding stations
visit_data_primary <- subset(visit_data_dual, visit_data_dual$perch == "primary")

#Add position x year identifier
visit_data_primary$position_year <- paste(visit_data_primary$position, visit_data_primary$year, sep = "_")

#Add info about previous visit at the same position
visit_data_primary$interval <- NULL

visit_data_primary <-visit_data_primary %>%
  arrange(position_year, start) %>%   
  group_by(position_year) %>%
  mutate(pre_JID = lag(JID),
         pre_start = lag(start),
         pre_end = lag(end),
         this_start = start,
         this_end = end,
         time_diff_pre = as.numeric(start - pre_end),
         pre_event_id = lag(event_id)) %>%
  ungroup()

overlaps <- visit_data_primary %>%
  arrange(position_year, start) %>% 
  group_by(position_year) %>%
  mutate(pre_end = lag(end),
         time_diff_pre = as.numeric(start - pre_end)) %>%
  filter(time_diff_pre < 0) %>%
  ungroup()

nrow(overlaps)
head(overlaps)

#Add info about the next visit at the same position
visit_data_primary <- visit_data_primary %>%
  arrange(position_year, start) %>%   
  group_by(position_year) %>%
  mutate(
    next_JID = lead(JID),
    next_start = lead(start),
    next_end = lead(end),
    time_diff_next = as.numeric(next_start - end),
    next_event_id = lead(event_id)
  ) %>%
  ungroup()

#Add column for self-follow: is the individual currently visiting the same individual that will visit the same location next
visit_data_primary$self_follow <- ifelse(visit_data_primary$JID == visit_data_primary$next_JID, "self", "not self")

table(visit_data_primary$self_follow)

#Subset data without self-follows
#visit_data_primary <- subset(visit_data_primary, !visit_data_primary$self_follow == "self")

#Parents of individual at primary perch
visit_data_primary$motherID <- LH_mother$MOTHER.ID[match(visit_data_primary$JID, LH_mother$ID)]
visit_data_primary$fatherID <- LH_father$FATHER.ID[match(visit_data_primary$JID, LH_father$ID)]

#Is the bird that visited just before the parent?
visit_data_primary$parent <- ifelse(visit_data_primary$pre_JID == visit_data_primary$motherID | visit_data_primary$pre_JID == visit_data_primary$fatherID, "parent", "no parent")
visit_data_primary$parent <- ifelse(is.na(visit_data_primary$parent), "no parent", visit_data_primary$parent)

#Age of the bird at primary perch
visit_data_primary$age_24 <- Ringed$min_age[match(visit_data_primary$JID, Ringed$ID)]
visit_data_primary$age <- ifelse(visit_data_primary$year == 2024, visit_data_primary$age_24, (visit_data_primary$age_24 -1))
table(visit_data_primary$age, visit_data_primary$parent)

#Jackdaw ID x position and JID x position x year combination 
visit_data_primary$JID_position <- paste(visit_data_primary$JID, visit_data_primary$position, sep = "_")
visit_data_primary$JID_position_year <- paste(visit_data_primary$JID, visit_data_primary$position, visit_data_primary$year, sep = "_")

#Subset primary perch visits by juveniles
visit_data_primary_juv <- subset(visit_data_primary, age == 1)

#First visit for each individual at each feeder 
visit_data_primary_first <- visit_data_primary %>%
  arrange(position, start) %>%   
  group_by(JID_position_year) %>%
  slice(1)

#Subset first feeder visits for juveniles
visit_data_primary_first_juv <- subset(visit_data_primary_first, age ==1)

#First visit just after parents: binary (0 = no, 1 = yes)
visit_data_primary_juv$parent <- ifelse(is.na(visit_data_primary_juv$parent), 0, ifelse(visit_data_primary_juv$parent == "parent", 1, 0))

#Subset instances where juveniles visit feeder for the first time just after parents
visit_data_primary_first_juv_parent <- subset(visit_data_primary_first_juv, parent == "parent")

#Signify "dual" juveniles that visited with their parents at least once 
visit_data_primary_first_juv_parent$dual_juv <- 1
visit_data_primary_first_juv$dual_juv <- visit_data_primary_first_juv_parent$dual_juv[match(visit_data_primary_first_juv$JID, visit_data_primary_first_juv_parent$JID)]
visit_data_primary_first_juv$dual_juv <- ifelse(is.na(visit_data_primary_first_juv$dual_juv), 0, 1)
table(visit_data_primary_first_juv$dual_juv)

#Signify specific instances where juveniles visited feeder for the first time and parent was there just before 
visit_data_primary_first_juv$parent <- ifelse(is.na(visit_data_primary_first_juv$parent), "no parent", visit_data_primary_first_juv$parent)
visit_data_primary_first_juv$parent <- ifelse(visit_data_primary_first_juv$parent == "parent", 1, 0)

#Subset "dual" juveniles that visited just after their parents at least once 
visit_data_primary_first_juv_dual <- subset(visit_data_primary_first_juv, dual_juv == 1)

#Subset "solo" juveniles that did not visit just after their parents
visit_data_primary_first_juv_solo <- subset(visit_data_primary_first_juv, !dual_juv == 1)

#Combine "dual"/"solo" juveniles with specific instances where parents visited just before; 0 = always solo, 1 = dual juvenile but solo at this feeder, 2 = dual at this feeder
visit_data_primary_first_juv$dual_juv_parent <- visit_data_primary_first_juv$dual_juv + visit_data_primary_first_juv$parent
visit_data_primary_first_juv$dual_juv_parent <- factor(visit_data_primary_first_juv$dual_juv_parent, levels = c(0, 1, 2))

#Add visit number at different feeders
visit_data_primary_visit_no <- visit_data_primary %>%
  count(JID_position)
visit_data_primary_first_juv$visit_no <- visit_data_primary_visit_no$n[match(visit_data_primary_first_juv$JID_position, visit_data_primary_visit_no$JID_position)]

#Add ringing date
visit_data_primary_first_juv$day_ringed <- LH_ring$day_ringed[match(visit_data_primary_first_juv$JID, LH_ring$ID)]
boxplot(visit_data_primary_first_juv$day_ringed ~ visit_data_primary_first_juv$dual_juv)
tapply(visit_data_primary_first_juv$day_ringed, visit_data_primary_first_juv$dual_juv, mean)

#Add nestbox
visit_data_primary_first_juv$box <- LH_ring$BOX[match(visit_data_primary_first_juv$JID, LH_ring$ID)]
visit_data_primary_first_juv$box[visit_data_primary_first_juv$box == ""] <- NA

#When did parents visit particular feeders for the first time?
visit_data_primary_first_juv$motherID_position_year <- paste(visit_data_primary_first_juv$motherID, visit_data_primary_first_juv$position, visit_data_primary_first_juv$year, sep = "_")
visit_data_primary_first_juv$fatherID_position_year <- paste(visit_data_primary_first_juv$fatherID, visit_data_primary_first_juv$position, visit_data_primary_first_juv$year, sep = "_")

visit_data_primary_first_juv$mother_position_first <- visit_data_primary_first$day[match(visit_data_primary_first_juv$motherID_position_year, visit_data_primary_first$JID_position_year)]
visit_data_primary_first_juv$father_position_first <- visit_data_primary_first$day[match(visit_data_primary_first_juv$fatherID_position_year, visit_data_primary_first$JID_position_year)]

#visit_data_primary_first_juv$mother_position_first[is.na(visit_data_primary_first_juv$mother_position_first)] <- 230 #if never visited, set after end of study period
#visit_data_primary_first_juv$father_position_first[is.na(visit_data_primary_first_juv$father_position_first)] <- 230 #if never visited, set after end of study period

visit_data_primary_first_juv$mother_days_before <- visit_data_primary_first_juv$day - visit_data_primary_first_juv$mother_position_first
visit_data_primary_first_juv$father_days_before <- visit_data_primary_first_juv$day - visit_data_primary_first_juv$father_position_first

visit_data_primary_first_juv$motherID_year <- paste(visit_data_primary_first_juv$motherID, visit_data_primary_first_juv$year, sep = "_")
visit_data_primary_first_juv$mother_visit_no_total <- visits$visit_number[match(visit_data_primary_first_juv$motherID_year, visits$JID_year)]

visit_data_primary_first_juv$fatherID_year <- paste(visit_data_primary_first_juv$fatherID, visit_data_primary_first_juv$year, sep = "_")
visit_data_primary_first_juv$father_visit_no_total <- visits$visit_number[match(visit_data_primary_first_juv$fatherID_year, visits$JID_year)]

#Add parent feeder use: 0 = never detected in that year, 1 = detected in that year but not at this position, 2 = detected this year and also at this position
visit_data_primary_first_juv$mother_feeder_use <- ifelse(is.na(visit_data_primary_first_juv$mother_visit_no_total) & is.na(visit_data_primary_first_juv$mother_days_before), 0, ifelse(!is.na(visit_data_primary_first_juv$mother_visit_no_total) & is.na(visit_data_primary_first_juv$mother_days_before), 1, 2))
visit_data_primary_first_juv$father_feeder_use <- ifelse(is.na(visit_data_primary_first_juv$father_visit_no_total) & is.na(visit_data_primary_first_juv$father_days_before), 0, ifelse(!is.na(visit_data_primary_first_juv$father_visit_no_total) & is.na(visit_data_primary_first_juv$father_days_before), 1, 2))
visit_data_primary_first_juv$parent_feeder_use <- visit_data_primary_first_juv$mother_feeder_use + visit_data_primary_first_juv$father_feeder_use

#Add age of individual just before 
visit_data_primary_first_juv$pre_age24 <- Ringed$min_age[match(visit_data_primary_first_juv$pre_JID, Ringed$ID)]
visit_data_primary_first_juv$pre_age <- ifelse(visit_data_primary_first_juv$year == 2024, visit_data_primary_first_juv$pre_age24, (visit_data_primary_first_juv$pre_age24 -1))

visit_data_primary_first$pre_JID_position_year <- paste(visit_data_primary_first$pre_JID, visit_data_primary_first$position_year, sep = "_")
visit_data_primary_first_juv$pre_JID_position_year <- paste(visit_data_primary_first_juv$pre_JID, visit_data_primary_first_juv$position_year, sep = "_")

visit_data_primary_first$pre_first_day <- visit_data_primary_first$day[match(visit_data_primary_first$pre_JID_position_year, visit_data_primary_first$JID_position_year)]
visit_data_primary_first_juv$pre_first_day <- visit_data_primary_first$pre_first_day[match(visit_data_primary_first_juv$pre_JID_position_year, visit_data_primary_first$pre_JID_position_year)]

#How many days before did the individual followed by the focal juvenile visit this position for the first time
visit_data_primary_first_juv$pre_days_before <- visit_data_primary_first_juv$day - visit_data_primary_first_juv$pre_first_day

#Juveniles following other juveniles 
visit_data_primary_first_juv$juv_juv <- ifelse(visit_data_primary_first_juv$age == 1 & visit_data_primary_first_juv$pre_age == 1, "juv_juv", "no_juv_juv")

#Subset of juveniles ringed before first juvenile started using feeders
visit_data_primary_first_juv2 <- subset(visit_data_primary_first_juv, day_ringed < 170)

#Just some quick summaries and visualisation
table(visit_data_primary_first_juv2$parent)
table(visit_data_primary_first_juv_dual$parent)
boxplot(visit_data_primary_first_juv_dual$day ~ visit_data_primary_first_juv_dual$parent)
table(visit_data_primary_first_juv_solo$parent)
mean(visit_data_primary_first_juv_solo$day)
table(visit_data_primary_first_juv$dual_juv_parent)

sum(visit_data_primary_first_juv2$mother_days_before > 0, na.rm = T)
sum(visit_data_primary_first_juv2$father_days_before > 0, na.rm = T)
sum(visit_data_primary_first_juv2$mother_days_before & visit_data_primary_first_juv2$father_days_before  < 0, na.rm = T)

table(visit_data_primary_first_juv2$parent_feeder_use)
table(visit_data_primary_first_juv2$mother_feeder_use)
table(visit_data_primary_first_juv2$father_feeder_use)

table(visit_data_primary_first_juv2$parent_feeder_use, visit_data_primary_first_juv2$dual_juv_parent)
table(visit_data_primary_first_juv2$parent_feeder_use[visit_data_primary_first_juv2$parent_feeder_use ==2], visit_data_primary_first_juv2$mother_feeder_use[visit_data_primary_first_juv2$parent_feeder_use== 2])

boxplot(visit_data_primary_first_juv2$day ~ visit_data_primary_first_juv2$juv_juv)
plot(visit_data_primary_first_juv2$pre_age, visit_data_primary_first_juv2$day)

#Add dual event identifier to total dataset ----

#Mark all dual events
dual_events$dual_event_key <- "dual_event"

#How many unique event ids at primary and secondary perches
dual_events_primary_event_ids <- dual_events  %>%  count(event_id)
dual_events_secondary_event_ids <- dual_events  %>%  count(secondary_event_id)

#Add dual events to full feeder dataset (dyadic feeding stations only)
visit_data_dual$dual_event_key1 <- dual_events$dual_event_key[match(visit_data_dual$event_id, dual_events$event_id)]
visit_data_dual$dual_event_key2 <- dual_events$dual_event_key[match(visit_data_dual$event_id, dual_events$secondary_event_id)]
visit_data_dual$dual_event_key <- ifelse(!is.na(visit_data_dual$dual_event_key1) | !is.na(visit_data_dual$dual_event_key2), "dual_event", "solo_event")

#Add dual events to primary perch dataset from section above
visit_data_primary$dual_event_key <- visit_data_dual$dual_event_key[match(visit_data_primary$event_id, visit_data_dual$event_id)]

#If no dual event, set event to "solo"
visit_data_primary$dual_event_key[is.na(visit_data_primary$dual_event_key)] <- "solo_event"

#Was event before a dual event
visit_data_primary$pre_dual_event_key <- visit_data_dual$dual_event_key[match(visit_data_primary$pre_event_id, visit_data_dual$event_id)]

#If no dual event, set relationship to "solo"
visit_data_primary$pre_dual_event_key[is.na(visit_data_primary$pre_dual_event_key)] <- "solo_event"

#Only keep visits at primary perches or at secondary perches if dual event
#visit_data_dual <- subset(visit_data_dual, perch == "primary" | dual_event_key == "dual_event")

#How many unique events
#visit_data_dual_event_ids <- visit_data_dual  %>%  count(event_id)
visit_data_primary_event_ids <- visit_data_primary  %>%  count(event_id)

#How many perch visits are part of dual events?
#table(visit_data_dual$dual_event_key)
table(visit_data_primary$dual_event_key)

#Merge all feeder visits at dyadic feeding stations with dual events dataset
#visit_data_dual2 <- merge(x = visit_data_dual, y = dual_events[, c("event_id", "arriv_diff", "relationship", "secondary_perch_id", "secondary_perch_JID", "dyad_ID")], by = "event_id", all = TRUE)
visit_data_primary2 <- merge(x = visit_data_primary, y = dual_events[, c("event_id", "secondary_event_id", "duration", "arriv_diff", "relationship", "secondary_perch_id", "secondary_perch_JID", "dyad_ID", "age_class_combination_simple")], by = "event_id", all = TRUE)

#If no dual event, set relationship to "solo"
visit_data_primary2$relationship[is.na(visit_data_primary2$relationship)] <- "solo"

#Add JID at primary feeder for those rows that were added by merging with dual event data
visit_data_primary2$JID <- ifelse(is.na(visit_data_primary2$JID), dual_events$primary_perch_JID[match(visit_data_primary2$event_id, dual_events$event_id)], visit_data_primary2$JID)

#Subset to only include visits at primary perches (visits at secondary perch are now logged due to merging with dual evens data)
#visit_data_dual2 <- subset(visit_data_dual2, perch == "primary")

#Finding N = 5540 dual events in full feeder visit data now
#table(visit_data_dual2$dual_event_key)
table(visit_data_primary2$dual_event_key)

#Adjust dual event key
#visit_data_primary2$dual_event_key <- ifelse(!is.na(visit_data_primary2$relationship), "dual_event", "solo_event")

#Add JID of queuer and relationship in previous event
visit_data_primary2 <- visit_data_primary2 %>%
  arrange(position_year, start) %>%   
  group_by(position_year) %>%
  mutate(
    pre_queuer_JID = lag(secondary_perch_JID),
    pre_relationship = lag(relationship),
    pre_dyad_ID = lag(dyad_ID),
    pre_age_class_combination_simple = lag(age_class_combination_simple),
    pre_duration = lag(duration),
    pre_arriv_diff = lag(arriv_diff)) %>%
  ungroup()

#Did the individual at the primary perch queue in the event before? Same if yes, different if someone else queued, solo if nobody queued before
visit_data_primary2$JID_pre_queuer_JID_same <- ifelse(visit_data_primary2$JID == visit_data_primary2$pre_queuer_JID, "same", "different") 
visit_data_primary2$JID_pre_queuer_JID_same <- ifelse(is.na(visit_data_primary2$JID_pre_queuer_JID_same), "solo", visit_data_primary2$JID_pre_queuer_JID_same)

#Did the individual at the primary perch queue in the event before? If there was a queueing individual in the previous event but the next visitor at primary perch is not the queuing individual - is it the one who was already at the primary perch or is it a third bird?
visit_data_primary2$JID_pre_queuer_JID_same_self_follow <- paste(visit_data_primary2$JID_pre_queuer_JID_same, visit_data_primary2$self_follow, sep = "_")

table(visit_data_primary2$JID_pre_queuer_JID_same)
table(visit_data_primary2$pre_relationship)
table(visit_data_primary2$JID_pre_queuer_JID_same, visit_data_primary2$pre_relationship)
table(visit_data_primary2$JID_pre_queuer_JID_same, visit_data_primary2$self_follow)
table(visit_data_primary2$pre_relationship, visit_data_primary2$self_follow, visit_data_primary2$JID_pre_queuer_JID_same)

#Add dual event key to juveniles' first visits: was the event before they first visited a feeder a dual even where they queued?
visit_data_primary_first_juv$dual_event_key <- visit_data_primary$pre_dual_event_key[match(visit_data_primary_first_juv$event_id, visit_data_primary$pre_event_id)]

length(unique(visit_data_primary_first_juv2$JID))
table(visit_data_primary_first_juv$dual_event_key, visit_data_primary_first_juv$dual_juv_parent)
table(visit_data_primary_first_juv$dual_event_key, visit_data_primary_first_juv$parent)
table(visit_data_primary_first_juv2$parent)

tapply(visit_data_primary_first_juv$time_diff_pre, visit_data_primary_first_juv$dual_juv_parent, mean, na.rm = T)
tapply(visit_data_primary_first_juv$time_diff_pre, visit_data_primary_first_juv$dual_juv_parent, sd, na.rm = T)
tapply(visit_data_primary_first_juv$time_diff_pre, visit_data_primary_first_juv$dual_juv_parent, median, na.rm = T)
tapply(visit_data_primary_first_juv$time_diff_pre, visit_data_primary_first_juv$parent, median, na.rm = T)

#Three contexts: bonded partner, other, solo
visit_data_primary2 <- visit_data_primary2 %>%
  mutate(relationship2 = case_when(
    relationship == "pair" ~ "bonded",
    relationship == "kin" ~ "bonded",
    relationship == "foraging associate" ~ "other",
    relationship == "neighbour" ~ "other",
    relationship == "other site" ~ "other",
    relationship == "site resident" ~ "other",
    relationship == "solo" ~ "solo"))

visit_data_primary2 <- visit_data_primary2 %>%
  mutate(pre_relationship2 = case_when(
    pre_relationship == "pair" ~ "bonded",
    pre_relationship == "kin" ~ "bonded",
    pre_relationship == "foraging associate" ~ "other",
    pre_relationship == "neighbour" ~ "other",
    pre_relationship == "other site" ~ "other",
    pre_relationship == "site resident" ~ "other",
    pre_relationship == "solo" ~ "solo")) 

visit_data_primary2$age_class_combination_simple <- ifelse(is.na(visit_data_primary2$age_class_combination_simple), "solo", visit_data_primary2$age_class_combination_simple)

#Displacement ----

visit_data_primary2$displaced <- ifelse(visit_data_primary2$time_diff_next < 2, 1, 0)

visit_data_primary2$JID_queuer_JID_displacer_same <- ifelse(visit_data_primary2$secondary_perch_JID == visit_data_primary2$next_JID, "same", "different")
visit_data_primary2$JID_queuer_JID_displacer_same <- ifelse(is.na(visit_data_primary2$JID_queuer_JID_displacer_same), "solo", visit_data_primary2$JID_queuer_JID_displacer_same)

table(visit_data_primary2$displaced, visit_data_primary2$JID_queuer_JID_displacer_same)
table(visit_data_primary2$displaced, visit_data_primary2$relationship)
table(visit_data_primary2$displaced, visit_data_primary2$relationship, visit_data_primary2$JID_queuer_JID_displacer_same)
table(visit_data_primary2$displaced, visit_data_primary2$dual_event_key)

#Only individuals that visited in three contexts ----
#(to compare contexts: with bonded adult, with another adult, solo)
visit_data_primary2$secondary_perch_age24 <- Ringed$age24[match(visit_data_primary2$secondary_perch_JID, Ringed$ID)]
visit_data_primary2$secondary_perch_age <- ifelse(visit_data_primary2$year == 2023, visit_data_primary2$secondary_perch_age24 - 1, visit_data_primary2$secondary_perch_age24)
visit_data_primary2$secondary_perch_age <- ifelse(is.na(visit_data_primary2$secondary_perch_age), "solo", visit_data_primary2$secondary_perch_age)

visit_data_primary2$next_age_24 <- Ringed$min_age[match(visit_data_primary2$next_JID, Ringed$ID)]
visit_data_primary2$next_age <- ifelse(visit_data_primary2$year == 2024, visit_data_primary2$next_age_24, (visit_data_primary2$next_age_24 -1))

#Remove cases where individual at secondary perch is juvenile 
visit_data_primary3 <- subset(visit_data_primary2, !visit_data_primary2$secondary_perch_age == 1)

visit_data_primary3 <- visit_data_primary3 %>%
  group_by(JID) %>%
  mutate(
    n_contexts = n_distinct(relationship2),
    all_contexts = n_contexts == n_distinct(visit_data_primary3$relationship2)
  ) %>%
  ungroup()

visit_data_primary3 <- subset(visit_data_primary3, visit_data_primary3$n_contexts == 3)

length(unique(visit_data_primary3$JID))

visit_data_primary4 <- subset(visit_data_primary3, !visit_data_primary3$JID_queuer_JID_displacer_same == "same")
visit_data_primary4 <- subset(visit_data_primary4, !visit_data_primary4$self_follow == "self")
visit_data_primary4 <- subset(visit_data_primary4, !visit_data_primary4$time_diff_next > 300)

visit_data_primary4$dyad_ID <- ifelse(!is.na(visit_data_primary4$dyad_ID), visit_data_primary4$dyad_ID, ifelse(visit_data_primary4$JID < visit_data_primary4$next_JID, paste(visit_data_primary4$JID, visit_data_primary4$next_JID, sep = "_"), paste(visit_data_primary4$next_JID, visit_data_primary4$JID, sep = "_"))) 

visit_data_primary4$next_age_24 <- Ringed$min_age[match(visit_data_primary4$next_JID, Ringed$ID)]
visit_data_primary4$next_age <- ifelse(visit_data_primary4$year == 2024, visit_data_primary4$next_age_24, (visit_data_primary4$next_age_24 -1))

visit_data_primary4$age_class <- ifelse(visit_data_primary4$age >1, "adult", "juvenile")
visit_data_primary4$next_age_class <- ifelse(visit_data_primary4$next_age >1, "adult", "juvenile")

visit_data_primary4$age_class_combination <- paste(visit_data_primary4$age_class, visit_data_primary4$next_age_class, sep = "_")


#Individual jackdaws data ----

#All individuals 
all_individuals <- as.data.frame(table(visit_data$JID))
all_individuals <- all_individuals %>% rename(JID = Var1)
all_individuals <- all_individuals %>% rename(visit_number = Freq)

all_individuals$sex <- LH_sex$SEX[match(all_individuals$JID, LH_sex$ID)]
all_individuals$sex <- ifelse(is.na(all_individuals$sex), "U", all_individuals$sex)

all_individuals$age23 <- Ringed$age23[match(all_individuals$JID, Ringed$ID)]
all_individuals$age24 <- Ringed$age24[match(all_individuals$JID, Ringed$ID)]


#Dual events individuals
dual_events_individuals <- as.data.frame(c(dual_events$initiator_JID, dual_events$joiner_JID))
dual_events_individuals <- dual_events_individuals %>% rename(JID = `c(dual_events$initiator_JID, dual_events$joiner_JID)`)

dual_events_individuals$sex <- LH_sex$SEX[match(dual_events_individuals$JID, LH_sex$ID)]
dual_events_individuals$sex <- ifelse(is.na(dual_events_individuals$sex), "U", dual_events_individuals$sex)

dual_events_individuals$age23 <- Ringed$age23[match(dual_events_individuals$JID, Ringed$ID)]
dual_events_individuals$age24 <- Ringed$age24[match(dual_events_individuals$JID, Ringed$ID)]

dual_events_individuals <- dual_events_individuals  %>%  
  add_count(JID, name = "dual_event_number")

dual_events_individuals <- dual_events_individuals %>%
  distinct()

table(dual_events_individuals$sex)

#Write datasets as CSV files ----
write.csv(dual_events,"dual_events.csv", row.names = FALSE)
write.csv(dual_events_dyads,"dual_events_dyads.csv", row.names = FALSE)
write.csv(visit_data_primary4,"visit_data_primary4.csv", row.names = FALSE)
write.csv(visit_data_primary_first_juv2,"visit_data_primary_first_juv2.csv", row.names = FALSE)
write.csv(visit_data_primary5,"visit_data_primary5.csv", row.names = FALSE)

#Read datasets for analyses ----
dual_events_dyads <- read.csv("dual_events_dyads.csv")
dual_events <- read.csv("dual_events.csv")
visit_data_primary4 <- read.csv("visit_data_primary4.csv")
visit_data_primary5 <- read.csv("visit_data_primary5.csv")
visit_data_primary_first_juv2 <- read.csv("visit_data_primary_first_juv2.csv") 
feedercoord <- read.csv("feedercoord.csv")


#(2) SUMMARY ----

#Getting an overview of the data, sample sizes etc. 

#Number of dual events
nrow(dual_events)

#Initiator at the primary perch
table(dual_events$initiator_perch)

#Duration of study period
max(visit_data$day[visit_data$year == 2023]) - min(visit_data$day[visit_data$year == 2023])
max(visit_data$day[visit_data$year == 2024]) - min(visit_data$day[visit_data$year == 2024])

#Number of individual jackdaws visiting any feeder
length(unique(visit_data$JID))

#Number of individual jackdaws engaging in dual event
length(unique(c(dual_events$initiator_JID, dual_events$joiner_JID)))

#Number of unique dyads engaging in dual event
length(unique(dual_events$dyad_ID))
length(unique(dual_events$dyad_ID[dual_events$day > 173]))

#Number of dyadic events during period where juveniles were visiting feeders
sum(dual_events$day > 173)

#How many dyads were pairs and kin
length(unique(dual_events$dyad_ID[dual_events$relationship == "pair"]))
length(unique(dual_events$dyad_ID[dual_events$relationship == "kin"]))

summary(dual_events)

#number of events per relationship type
table(dual_events$relationship)

#number of unique dyads where by AB and BA are considered the same
length(table(dual_events$dyad_ID)) #1001 unique dyads

#number of dyads where by AB and BA are considered separately
length(table(dual_events$initiator_joiner_ID)) #1283 unique dyad combinations

#number of different initiatiors and joiners
length(table(dual_events$initiator_JID)) #192
length(table(dual_events$joiner_JID)) #246

length(table(dual_events_pair$initiator)) #43
length(table(dual_events_pair$joiner)) #43

table(dual_events$initiator)
table(dual_events$joiner)

table(dual_events_individuals$sex) #52 F, 74 M, 117 NA

table(dual_events$initiator_sex) #676 F, 1531 M, 683 NA
table(dual_events$joiner_sex) #764 F, 1337 M, 789 NA
binom.test(676, 2207, 0.5)
binom.test(764, 2101, 0.5)

#Subset sex combination 
dual_events_sex <- subset(dual_events[, c("initiator_sex", "joiner_sex")])
dual_events_sex <- dual_events_sex %>%
  mutate(initiator_sex =case_when(
    is.na(initiator_sex) ~ "N",
    initiator_sex == "F" ~ "F",
    initiator_sex == "M" ~ "M"))
dual_events_sex <- dual_events_sex %>%
  mutate(joiner_sex =case_when(
    is.na(joiner_sex) ~ "N",
    joiner_sex == "F" ~ "F",
    joiner_sex == "M" ~ "M"))

table(dual_events$sex_combination)

#Chi Squared test about sex combinations
sex_comb_matrix <- as.table(rbind(c(77, 355, 127), c(456, 516, 279), c(91, 188, 135)))
dimnames(sex_comb_matrix) <- list(initiator_sex=c("Female","Male","Nonsexed"),joiner_sex=c("Female","Male","Nonsexed"))
sex_comb_matrix
chisq.test(sex_comb_matrix)
chisq.posthoc.test(sex_comb_matrix, method = "bonferroni")

table(dual_events$initiator_feeder)
table(dual_events$joiner_feeder)

#Sex combinations for pairs
table(dual_events$sex_combination)
table(dual_events$pair) #1651 non-pair, 429 pair
table(dual_events_pair$sex_combination) #pairs: 143 FM, 256 MF
table(dual_events_pair$initiator_sex)

#Kin combinations
table(dual_events$kin) #221 kin, 699 non-kin, 1898 NA
table(dual_events$parent_offspring) #133 parent offspring
table(dual_events$offspring_parent) #53 offspring parent
table(dual_events$kin, dual_events$age_class_combination)
table(dual_events$relationship) #186 kin, 429 pair, 2275 other

#Initiator and joiner duration for pairs per sex
aggregate(dual_events_pair_primary$initiator_duration, by=list(dual_events_pair_primary$sex_combination), FUN=mean)
aggregate(dual_events_pair_primary$joiner_duration, by=list(dual_events_pair_primary$sex_combination), FUN=mean)

#Dual event/overlap duration per relationship type
aggregate(dual_events$overlap_duration, by=list(dual_events$relationship), FUN=mean) #kin = 10.68, other = 3.49, pair = 11.69
aggregate(dual_events$overlap_duration, by=list(dual_events$relationship), FUN=sd) #kin = 11.62, other = 5.47, pair = 13.36

aggregate(dual_events$overlap_duration, by=list(dual_events$kin_age), FUN=mean) 
aggregate(dual_events$overlap_duration, by=list(dual_events$kin_age), FUN=sd) 

#Initiator and joiner duration per relationship type
aggregate(dual_events_primary$initiator_duration, by=list(dual_events_primary$relationship), FUN=mean) # kin = 35.53, other = 37.04, pair = 40.99
aggregate(dual_events_primary$initiator_duration, by=list(dual_events_primary$relationship), FUN=sd) # kin = 35.53, other = 37.04, pair = 40.99

aggregate(dual_events_primary$joiner_duration, by=list(dual_events_primary$relationship), FUN=mean) #kin = 16.86, other = 4.21, pair = 16.66
aggregate(dual_events_primary$joiner_duration, by=list(dual_events_primary$relationship), FUN=sd) #kin = 16.86, other = 4.21, pair = 16.66

#Initiator and joiner duration per sex combination
aggregate(dual_events_primary$initiator_duration, by=list(dual_events_primary$sex_combination), FUN=mean)
aggregate(dual_events_primary$joiner_duration, by=list(dual_events_primary$sex_combination), FUN=mean)
aggregate(dual_events_primary$initiator_duration, by=list(dual_events_primary$sex_combination_simple), FUN=mean)
aggregate(dual_events_primary$joiner_duration, by=list(dual_events_primary$sex_combination_simple), FUN=mean)

#Dual event duration per pair ID
aggregate(dual_events$duration, by=list(dual_events$pair_ID), FUN=mean)

#Average initiator duration
mean(dual_events$initiator_duration)
sd(dual_events$initiator_duration)

#Average overlap duration
mean(dual_events$overlap_duration)
sd(dual_events$overlap_duration)
median(dual_events$overlap_duration)

mean(dual_events$prop_dual_overlap)
sd(dual_events$prop_dual_overlap)

boxplot(dual_events$overlap_duration ~ dual_events$site)
boxplot(dual_events$arriv_diff ~ dual_events$site)
boxplot(dual_events_primary$initiator_duration ~ dual_events_primary$site)

boxplot(visit_data$visit_duration ~ visit_data$site)


#(3) STATISTICAL ANALYSIS ----

#1. Social preference ----

tapply(dual_events_dyads$visit_number_dyad, dual_events_dyads$relationship, mean)

dual_events_dyads$obs <- row.names(dual_events_dyads)

#Bayesian model 

#Model just with prior
default_prior()

social_bias_brm1_prior <- brm(visit_number_dyad | trunc(lb = 1) ~ relationship +
                           age_class_combination_simple + (1|position) +
                           (1|mm(initiator_JID, joiner_JID)), 
                           data = dual_events_dyads, 
                           family = poisson(link = "log"),
                           prior = c(
                             prior(normal(log(2.5), 1), class = "Intercept"),  
                             prior(normal(0, 0.5), class = "b"),
                             prior(exponential(1), class = "sd")),              
                           sample_prior = "only")

#Prior predictive checks 
pp_check(social_bias_brm1_prior, ndraws = 100) +
  xlim(0, 100)

pp_check(social_bias_brm1_prior, ndraws = 100) +
  scale_x_log10()

yrep <- posterior_predict(social_bias_brm1_prior)
quantile(yrep, 0.90)  # 99.9th percentile
max(yrep)              # largest simulated count

#Model

#Entire study period

#Poisson
social_bias_brm1_poisson <- brm(visit_number_dyad | trunc(lb = 1) ~ relationship +
                                age_class_combination_simple + (1|position) +
                                (1|mm(initiator_JID, joiner_JID)), 
                              data = dual_events_dyads, 
                              family = poisson(link = "log"),
                              prior = c(
                                prior(normal(log(2.5), 1), class = "Intercept"),  
                                prior(normal(0, 0.5), class = "b"),
                                prior(exponential(1), class = "sd")),
                              control = list(adapt_delta = 0.99, max_treedepth = 15),
                              save_pars = save_pars(all = TRUE)
)
   
#Negative binomial 
social_bias_brm1_negbinomial <- brm(visit_number_dyad | trunc(lb = 1) ~ relationship +
                          age_class_combination_simple + (1|position) +
                          (1|mm(initiator_JID, joiner_JID)), 
                        data = dual_events_dyads, 
                        family = negbinomial(link = "log"),
                        prior = c(
                          prior(normal(log(2.5), 1), class = "Intercept"),  
                          prior(normal(0, 0.5), class = "b"),
                          prior(exponential(1), class = "sd"),
                          prior(gamma(2, 0.1), class = "shape")),
)        

#Poisson with observation-level random effect
social_bias_brm1_poisson_obs <- brm(visit_number_dyad | trunc(lb = 1) ~ relationship +
                                  age_class_combination_simple + (1|position) + (1|obs) +
                                  (1|mm(initiator_JID, joiner_JID)), 
                                data = dual_events_dyads, 
                                family = poisson(link = "log"),
                                prior = c(
                                  prior(normal(log(2.5), 1), class = "Intercept"),  
                                  prior(normal(0, 0.5), class = "b"),
                                  prior(exponential(1), class = "sd")),
                                control = list(adapt_delta = 0.99, max_treedepth = 15),
                                save_pars = save_pars(all = TRUE)
)


#iter = 6000, warmup = 2000,
#control = list(adapt_delta = 0.99))

#Only period in which juveniles used feeders 
#Entire study period
social_bias_brm1_juv_poisson <- brm(visit_number_dyad | trunc(lb = 1) ~ relationship +
                                  age_class_combination_simple + (1|position) +
                                  (1|mm(initiator_JID, joiner_JID)), 
                                data = dual_events_dyads_juv, 
                                family = poisson(link = "log"),
                                prior = c(
                                  prior(normal(log(2.5), 1), class = "Intercept"),  
                                  prior(normal(0, 0.5), class = "b"),
                                  prior(exponential(1), class = "sd")),
)

social_bias_brm1 <- social_bias_brm1_poisson
social_bias_brm1 <- social_bias_brm1_negbinomial
social_bias_brm1 <- social_bias_brm1_poisson_obs

#Model summary
summary(social_bias_brm1_poisson)
summary(social_bias_brm1_poisson_obs)

#Check collinearity
check_collinearity(social_bias_brm1)

#Posterior predictive checks
pp_check(social_bias_brm1, type = "dens_overlay", ndraws = 100) 
pp_check(social_bias_brm1_poisson, type = "dens_overlay", ndraws = 100) + scale_x_log10()
pp_check(social_bias_brm1_poisson_obs, type = "dens_overlay", ndraws = 100) + scale_x_log10()

pp_check(social_bias_brm1, type = "dens_overlay", ndraws = 100) +
  scale_x_log10()

pp_check(social_bias_brm1, type = "stat", stat = "var", ndraws = 100)
pp_check(social_bias_brm1, type = "bars", ndraws = 100)
pp_check(social_bias_brm1, type = "hist", ndraws = 100)

y <- model.frame(social_bias_brm1)[[1]]
yrep <- posterior_predict(social_bias_brm1)

#Plot model
plot(social_bias_brm1)
launch_shinystan(social_bias_brm1)

#Posterior distribution
as_draws_df(social_bias_brm1)
mcmc_areas(social_bias_brm1)
mcmc_intervals(social_bias_brm1)

#Evaluation and interpretation
loo(social_bias_brm1, moment_match = TRUE)
loo(social_bias_brm1_poisson, social_bias_brm1_poisson_obs)
conditional_effects(social_bias_brm1)
conditional_effects(social_bias_brm1,effects = "relationship")
bayes_R2(social_bias_brm1)

#Sensitivity analysis and prior checks 
prior_summary(social_bias_brm1)

#Extract predictions
fitted(social_bias_brm1, scale = "response")

#Calculate raw mean, median, and SD per group 
tapply(dual_events_dyads$visit_number_dyad, dual_events_dyads$relationship, FUN = mean)
tapply(dual_events_dyads$visit_number_dyad, dual_events_dyads$relationship, FUN = median)
tapply(dual_events_dyads$visit_number_dyad, dual_events_dyads$relationship, FUN = sd)

#Posterior marginal 
dual_events_dyads$relationship <- factor(dual_events_dyads$relationship, levels = c("other site" , "site resident", "foraging associate", "neighbour", "kin", "pair"))

posterior_marginal <- dual_events_dyads %>%
  add_epred_draws(
    object = social_bias_brm1,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_relationship <- posterior_marginal %>%
  group_by(.draw, relationship) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_summary <- posterior_draw_relationship %>%
  group_by(relationship) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_relationship %>%
  compare_levels(epred, by = relationship)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = relationship,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

#Pairwise contrasts for Supp Mat using emmeans
contrasts <- social_bias_brm1 %>%
  emmeans(~ relationship) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  median_qi()
contrasts

#Hypothesis tests 
hypothesis(social_bias_brm1, "relationshippair > 0")
hypothesis(social_bias_brm1, "relationshippair - relationshipkin > 0")
hypothesis(social_bias_brm1, "relationshippair - relationshipneighbour > 0")
hypothesis(social_bias_brm1, "relationshippair - relationshipsite resident > 0")
hypothesis(social_bias_brm1, "relationshippair - relationshipother site > 0")

hypothesis(social_bias_brm1, "relationshipkin > 0")
hypothesis(social_bias_brm1, "relationshipkin - relationshipneighbour > 0")
hypothesis(social_bias_brm1, "relationshipkin - relationshipsite resident > 0")
hypothesis(social_bias_brm1, "relationshipkin - relationshipother site > 0")

hypothesis(social_bias_brm1, "relationshipneighbour > 0")
hypothesis(social_bias_brm1, "relationshipneighbour - relationshipsite resident > 0")
hypothesis(social_bias_brm1, "relationshipneighbour - relationshipother site > 0")

hypothesis(social_bias_brm1, "relationshipsite resident < 0")
hypothesis(social_bias_brm1, "relationshipother site < 0")

hypothesis(social_bias_brm1, "relationshipsite resident - relationshipother site > 0")

#Halfeye plot
social_bias_halfeye_plot <- ggplot(posterior_draw_relationship, aes(y = relationship, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "grey") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")
  ) +
  labs(x = "Dyadic events per dyad", y = "Relationship") +
  xlim(0, 10)
social_bias_halfeye_plot

social_bias_halfeye_plot_not_pairs <- ggplot(posterior_draw_relationship, aes(y = relationship, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "grey") +
  theme_few(base_size = 18) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")
  ) +
  labs(x = "Dyadic events per dyad", y = "Relationship") +
  xlim(2, 2.5)
social_bias_halfeye_plot_not_pairs

posterior_draw_relationship_pairs <- subset(posterior_draw_relationship, posterior_draw_relationship$relationship == "pair")

social_bias_halfeye_plot_pairs <- ggplot(posterior_draw_relationship_pairs, aes(y = relationship, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "grey") +
  theme_few(base_size = 40) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")
  ) +
  labs(x = "Dyadic events per dyad", y = "Relationship") +
  xlim(2, 9)
social_bias_halfeye_plot

#2. Coordination ----

#Bayesian model 

#Model just with prior
default_prior()

options(contrasts = c("contr.sum", "contr.poly"))

coordination_brm1_prior <- brm(arriv_diff ~ relationship +
                                age_class_combination_simple + (1|position) +
                                (1|mm(initiator_JID, joiner_JID)) + (1|dyad_ID), 
                              data = dual_events, 
                              family = Gamma(link = "log"),
                              prior = c(
                                  prior(normal(log(25), 0.25), class = "Intercept", ub = log(500)),
                                  prior(normal(0, 0.5), class = "b"),
                                  prior(exponential(1.5), class = "sd"),
                                  prior(gamma(2, 0.5), class = "shape")),
                              iter = 6000, warmup = 2000,
                              control = list(adapt_delta = 0.99, max_treedepth = 15),
                              sample_prior = "only")

                            

#Prior predictive checks 
pp_check(coordination_brm1_prior, ndraws = 100) +
  xlim(0, 50)

pp_check(coordination_brm1_prior, type="dens_overlay", ndraws = 100) + 
  scale_x_log10()

pp_check(coordination_brm1_prior, type = "ecdf_overlay", ndraws = 100) +
  xlim(0, 100)

yrep <- posterior_predict(coordination_brm1_prior)
quantile(yrep, 0.99)  # 99.9th percentile
max(yrep)              # largest simulated count

yrep <- posterior_predict(coordination_brm1_prior, draws = 100)


#Model
coordination_brm1<- brm(arriv_diff ~ relationship +
                                 age_class_combination_simple + (1|position) +
                                 (1|mm(initiator_JID, joiner_JID)) + (1|dyad_ID), 
                               data = dual_events, 
                               family = Gamma(link = "log"),
                               prior = c(
                                 prior(normal(log(25), 0.25), class = "Intercept", ub = log(500)),
                                 prior(normal(0, 0.5), class = "b"),
                                 prior(exponential(1.5), class = "sd"),
                                 prior(gamma(2, 0.5), class = "shape"))
                        
)
                               #iter = 6000, warmup = 2000,
                               #control = list(adapt_delta = 0.99, max_treedepth = 15))

#Model summary
summary(coordination_brm1)

#Check collinearity
check_collinearity(coordination_brm1)

#Posterior predictive checks
pp_check(coordination_brm1, type = "dens_overlay", ndraws = 100) +
  xlim(0, 200)

pp_check(coordination_brm1, type = "dens_overlay", ndraws = 100) +
  scale_x_log10()
  
pp_check(coordination_brm1, type = "dens_overlay", ndraws = 100) 
  
pp_check(coordination_brm1, type = "stat", stat = "var", ndraws = 1000)
pp_check(coordination_brm1, type = "bars", ndraws = 1000)
pp_check(coordination_brm1, type = "hist", ndraws = 1000)

y <- model.frame(coordination_brm1)[[1]]
yrep <- posterior_predict(coordination_brm1)

# Mean variance ratio across posterior draws:
dispersion <- apply(yrep, 1, function(x) var(x - y))
mean(dispersion)

#Plot model
plot(coordination_brm1)
launch_shinystan(coordination_brm1)

#Posterior distribution
as_draws_df(coordination_brm1)
mcmc_areas(coordination_brm1)
mcmc_intervals(coordination_brm1)

#First e.valuation and interpretation
loo(coordination_brm1, moment_match = TRUE)
fitted(coordination_brm1, scale = "response")
conditional_effects(coordination_brm1)
conditional_effects(coordination_brm1,effects = "relationship")
bayes_R2(coordination_brm1)
emmeans(coordination_brm1, ~ relationship, type = "response")
emmeans(coordination_brm1, ~ relationship, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(coordination_brm1)

#Extract predictions
fitted(coordination_brm1)

#Posterior marginal 
dual_events$relationship <- factor(dual_events$relationship, levels = c("other site" , "site resident", "foraging associate", "neighbour", "kin", "pair"))

posterior_marginal <- dual_events %>%
  add_epred_draws(
    object = coordination_brm1,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_relationship <- posterior_marginal %>%
  group_by(.draw, relationship) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_summary <- posterior_draw_relationship %>%
  group_by(relationship) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_relationship %>%
  compare_levels(epred, by = relationship)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = relationship,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

#Pairwise contrasts for Supp Mat using emmeans
contrasts <- coordination_brm1 %>%
  emmeans(~ relationship) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  median_qi()
contrasts

#Hypothesis tests 
hypothesis(coordination_brm1, "relationshippair > 0")
hypothesis(coordination_brm1, "relationshippair - relationshipkin > 0")
hypothesis(coordination_brm1, "relationshippair - relationshipneighbour > 0")
hypothesis(coordination_brm1, "relationshippair - relationshipsite resident > 0")
hypothesis(coordination_brm1, "relationshippair - relationshipother site > 0")

hypothesis(coordination_brm1, "relationshipkin > 0")
hypothesis(coordination_brm1, "relationshipkin - relationshipneighbour > 0")
hypothesis(coordination_brm1, "relationshipkin - relationshipsite resident > 0")
hypothesis(coordination_brm1, "relationshipkin - relationshipother site > 0")

hypothesis(coordination_brm1, "relationshipneighbour > 0")
hypothesis(coordination_brm1, "relationshipneighbour - relationshipsite resident > 0")
hypothesis(coordination_brm1, "relationshipneighbour - relationshipother site > 0")

hypothesis(coordination_brm1, "relationshipsite resident < 0")
hypothesis(coordination_brm1, "relationshipother site < 0")

hypothesis(coordination_brm1, "relationshipsite resident - relationshipother site > 0")

#Violin plot of raw data
dual_events$relationship <- factor(dual_events$relationship, levels = c("pair", "kin", "neighbour", "foraging associate", "site resident", "other site"))

coordination_violin_plot <- ggplot(dual_events, aes(x = relationship, y = arriv_diff)) +
  geom_violin(alpha = 1, position = position_dodge(width = 0.8), width = 0.7, fill = "white") +
  geom_jitter(
    aes(x = relationship, y = arriv_diff), 
    data = dual_events, 
    color = "black",
    size = 2,
    alpha = 0.05,
    inherit.aes = FALSE,
    width = 0.2) +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(size = 14, family = "Garamond")) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.8), color = "white") +
  labs(x = "Relationship", y = "Latency of arrival (s)") +
  scale_x_discrete(guide = guide_axis(angle = 15))
coordination_violin_plot

#Halfeye plot of model predictions
coordination_halfeye_plot <- ggplot(posterior_draw_relationship, aes(y = relationship, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "grey") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Latency of arrival (s)", y = "Relationship") +
  xlim(c(0, 50))
coordination_halfeye_plot

#Panel 
coordination_plot <- ggarrange(coordination_violin_plot, coordination_halfeye_plot, ncol = 2, nrow = 1,  widths = c(1, 1))
coordination_plot

#3. Social tolerance and support ----

#3.1. Social tolerance ----

#Bayesian model 

#Model just with prior
default_prior()

options(contrasts = c("contr.sum", "contr.poly"))

#Censoring durations detected as 0 

#Remove constant added above
#dual_events$duration <- dual_events$duration - 0.5

# Make a new column for modelling
dual_events$duration2 <- dual_events$duration

# Set censored values to an approximate detection limit
dual_events$duration2[dual_events$duration == 0] <- 1
summary(dual_events$duration2)

# Create a censoring indicator
dual_events$duration_cens <- ifelse(dual_events$duration == 0, "left", "none")
table(dual_events$duration_cens)

social_tolerance_brm1_prior <- brm(duration2 | cens(duration_cens) ~ relationship +
                                 age_class_combination_simple + (1|position) +
                                 (1|mm(initiator_JID, joiner_JID)), 
                               data = dual_events, 
                               family = Gamma(link = "log"),
                               prior = c(
                                 prior(normal(log(6), 0.5), class = "Intercept"),  
                                 prior(normal(0, 0.3), class = "b"),
                                 prior(exponential(1), class = "sd"),
                                 prior(gamma(2, 0.5), class = "shape")),
                               iter = 6000, warmup = 2000,
                               control = list(adapt_delta = 0.99, max_treedepth = 15),
                               sample_prior = "only")

#Prior predictive checks 
pp_check(social_tolerance_brm1_prior, ndraws = 100) +
  xlim(0, 20)

pp_check(social_tolerance_brm1_prior, type="dens_overlay", ndraws = 100) + 
  scale_x_log10()

yrep <- posterior_predict(social_tolerance_brm1_prior)
quantile(yrep, 0.85)  # 99.9th percentile
max(yrep)              # largest simulated count

#Model
dual_events$arriv_diff_z <- scale(dual_events$arriv_diff)

social_tolerance_brm1_gamma <- brm(duration2 | cens(duration_cens) ~ relationship + arriv_diff_z + 
                              age_class_combination_simple + (1|position) + (1|dyad_ID) +
                              (1|mm(initiator_JID, joiner_JID)), 
                              data = dual_events,
                              family = Gamma(link = "log"),
                              prior = c(
                                prior(normal(log(6), 0.5), class = "Intercept"),
                                prior(normal(0, 0.3), class = "b"),
                                prior(exponential(1), class = "sd"),
                                prior(gamma(2, 0.5), class = "shape")
                              ),
                              control = list(adapt_delta = 0.999, max_treedepth = 15),
                              save_pars = save_pars(all = TRUE)
)

                              #iter = 6000, warmup = 2000,
                              #control = list(adapt_delta = 0.99, max_treedepth = 15),


social_tolerance_brm1_lognormal <- brm(duration2 | cens(duration_cens) ~ relationship + arriv_diff +
                              age_class_combination_simple + (1|position) + (1|dyad_ID) +
                              (1|mm(initiator_JID, joiner_JID)),
                              data = dual_events,
                              family = lognormal(link = "identity"),
                              prior = c(
                                prior(normal(log(6), 0.5), class = "Intercept"),
                                prior(normal(0, 0.3), class = "b"),
                                prior(exponential(1), class = "sd"),
                                prior(exponential(1), class = "sigma")
                              ),
                              iter = 4000, warmup = 2000,
                              control = list(adapt_delta = 0.95, max_treedepth = 12),
                              save_pars = save_pars(all = TRUE)
)

social_tolerance_brm1 <- social_tolerance_brm1_gamma

#Model summary
summary(social_tolerance_brm1)
summary(social_tolerance_brm1_gamma)
summary(social_tolerance_brm1_lognormal)

#Check collinearity
check_collinearity(social_tolerance_brm1)

social_tolerance_brm1_gamma_loo <- loo(social_tolerance_brm1_gamma, moment_match = TRUE)
social_tolerance_brm1_lognormal_loo <- loo(social_tolerance_brm1_lognormal, moment_match = TRUE)
loo_compare(social_tolerance_brm1_gamma_loo, social_tolerance_brm1_lognormal_loo)

#Posterior predictive checks
pp_check(social_tolerance_brm1, type = "dens_overlay", ndraws = 100) 
  
pp_check(social_tolerance_brm1, type = "dens_overlay", ndraws = 100) +
  scale_x_log10()

pp_check(social_tolerance_brm1, type = "dens_overlay", ndraws = 100) +
  xlim(0, 50)

pp_check(social_tolerance_brm1_lognormal, type = "dens_overlay", ndraws = 100) +
  xlim(0, 100)

pp_check(social_tolerance_brm1, type = "dens_overlay", obs_args = list(col = "black"))

pp_check(social_tolerance_brm1, type = "stat", stat = "var", ndraws = 1000)
pp_check(social_tolerance_brm1, type = "bars", ndraws = 1000)
pp_check(social_tolerance_brm1, type = "hist", ndraws = 1000)

y <- model.frame(social_tolerance_brm1)[[1]]
yrep <- posterior_predict(social_tolerance_brm1)
# Mean variance ratio across posterior draws:
dispersion <- apply(yrep, 1, function(x) var(x - y))
mean(dispersion)

#Plot model
plot(social_tolerance_brm1)
launch_shinystan(social_tolerance_brm1)

#Posterior distribution
as_draws_df(social_tolerance_brm1)
mcmc_areas(social_tolerance_brm1)
mcmc_intervals(social_tolerance_brm1)

#First evaluation and interpretation
loo(social_tolerance_brm1, moment_match = TRUE)
fitted(social_tolerance_brm1, scale = "response")
conditional_effects(social_tolerance_brm1)
conditional_effects(social_tolerance_brm1,effects = "relationship")
bayes_R2(social_tolerance_brm1)
emmeans(social_tolerance_brm1, ~ relationship, type = "response")
emmeans(social_tolerance_brm1, ~ relationship, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(social_tolerance_brm1)

#Extract predictions
fitted(social_tolerance_brm1)

#Posterior marginal 
dual_events$relationship <- factor(dual_events$relationship, levels = c("other site" , "site resident", "foraging associate", "neighbour", "kin", "pair"))

posterior_marginal <- dual_events %>%
  add_epred_draws(
    object = social_tolerance_brm1,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_relationship <- posterior_marginal %>%
  group_by(.draw, relationship) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_summary <- posterior_draw_relationship %>%
  group_by(relationship) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_relationship %>%
  compare_levels(epred, by = relationship)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = relationship,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

#Pairwise contrasts for Supp Mat using emmeans
contrasts <- coordination_brm1 %>%
  emmeans(~ relationship) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  median_qi()
contrasts

hypothesis(social_tolerance_brm1, "relationshippair > 0")
hypothesis(social_tolerance_brm1, "relationshippair - relationshipkin > 0")
hypothesis(social_tolerance_brm1, "relationshippair - relationshipneighbour > 0")
hypothesis(social_tolerance_brm1, "relationshippair - relationshipsite resident > 0")
hypothesis(social_tolerance_brm1, "relationshippair - relationshipother site > 0")

hypothesis(social_tolerance_brm1, "relationshipkin > 0")
hypothesis(social_tolerance_brm1, "relationshipkin - relationshipneighbour > 0")
hypothesis(social_tolerance_brm1, "relationshipkin - relationshipsite resident > 0")
hypothesis(social_tolerance_brm1, "relationshipkin - relationshipother site > 0")

hypothesis(social_tolerance_brm1, "relationshipneighbour > 0")
hypothesis(social_tolerance_brm1, "relationshipneighbour - relationshipsite resident > 0")
hypothesis(social_tolerance_brm1, "relationshipneighbour - relationshipother site > 0")

hypothesis(social_tolerance_brm1, "relationshipsite resident > 0")
hypothesis(social_tolerance_brm1, "relationshipother site > 0")

hypothesis(social_tolerance_brm1, "relationshipsite resident - relationshipother site > 0")

#Violin plot of raw data
dual_events$relationship <- factor(dual_events$relationship, levels = c("pair", "kin", "neighbour", "foraging associate", "site resident", "other site"))

social_tolerance_violin_plot <- ggplot(dual_events, aes(x = relationship, y = duration)) +
  geom_violin(alpha = 1, position = position_dodge(width = 0.8), width = 0.7, fill = "white") +
  geom_jitter(
    aes(x = relationship, y = duration), 
    data = dual_events, 
    color = "black",
    size = 2,
    alpha = 0.05,
    inherit.aes = FALSE,
    width = 0.2) +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(size = 14, family = "Garamond")) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.8), color = "white") +
  labs(x = "Relationship", y = "Duration of dyadic event (s)") +
  scale_x_discrete(guide = guide_axis(angle = 15))
social_tolerance_violin_plot

#Halfeye plot of model predictions
social_tolerance_halfeye_plot <- ggplot(posterior_draw_relationship, aes(y = relationship, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "grey") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Duration of dyadic event (s)", y = "Relationship") +
  xlim(c(0, 20))
social_tolerance_halfeye_plot

#Panel 
social_tolerance_plot <- ggarrange(social_tolerance_violin_plot, social_tolerance_halfeye_plot, ncol = 2, nrow = 1,  widths = c(1, 1))
social_tolerance_plot

#Model prediction summary table 
pred_summary <- posterior_summary %>%
  group_by(relationship) %>%
  summarise(
    median_pred = median(.epred),
    lower_95 = quantile(.epred, 0.025),
    upper_95 = quantile(.epred, 0.975),
    lower_50 = quantile(.epred, 0.25),
    upper_50 = quantile(.epred, 0.75)
  )
pred_summary


#3.2. Forager benefit during dyadic events ----

#Bayesian model 

#Model just with prior
default_prior()

#Censoring durations detected as 0 

#Remove constant added above
#dual_events$duration <- dual_events$duration - 0.5

# Make a new column for modeling
dual_events$primary_duration2 <- dual_events$primary_duration

# Set censored values to an approximate detection limit
dual_events$primary_duration2[dual_events$primary_duration == 0] <- 1
summary(dual_events$primary_duration2)

# Create a censoring indicator
dual_events$primary_duration_cens <- ifelse(dual_events$primary_duration == 0, "left", "none")
table(dual_events$primary_duration_cens)

#Scale latency of arrival
dual_events$arriv_diff_z <- scale(dual_events$arriv_diff)

#Subset data: only include events where the initiator is at the primary perch
dual_events2 <- subset(dual_events, initiator_perch == "primary")

#Full dataset, censored response
social_tolerance_brm2_prior <- brm(primary_duration2 | cens(primary_duration_cens) ~ relationship +
                                     age_class_combination_simple + arriv_diff_z + (1|position) +
                                     (1|mm(initiator_JID, joiner_JID)), 
                                   data = dual_events, 
                                   family = Gamma(link = "log"),
                                   prior = c(
                                     prior(normal(log(35), 0.5), class = "Intercept"),
                                     prior(normal(0, 0.3), class = "b"),
                                     prior(exponential(1), class = "sd"),
                                     prior(gamma(2, 0.5), class = "shape")),
                                   sample_prior = "only")

#Subset of data where initiator at primary perch, no censored response
social_tolerance_brm2_prior <- brm(primary_duration ~ relationship +
                                     age_class_combination_simple + arriv_diff_z + (1|position) +
                                     (1|mm(initiator_JID, joiner_JID)), 
                                   data = dual_events2, 
                                   family = Gamma(link = "log"),
                                   prior = c(
                                     prior(normal(log(35), 0.5), class = "Intercept"),
                                     prior(normal(0, 0.3), class = "b"),
                                     prior(exponential(1), class = "sd"),
                                     prior(gamma(2, 0.5), class = "shape")),
                                   sample_prior = "only")

#Prior predictive checks 
pp_check(social_tolerance_brm2_prior, ndraws = 100) +
  xlim(0.5, 100)

pp_check(social_tolerance_brm2_prior, type="dens_overlay", ndraws = 100) + 
  scale_x_log10()

yrep <- posterior_predict(social_tolerance_brm2_prior)
quantile(yrep, 0.85)  # 99.9th percentile
max(yrep)              # largest simulated count

#Model

#Full dataset, censored response
social_tolerance_brm2_gamma <- brm(primary_duration2 | cens(primary_duration_cens) ~ 
                                   relationship + age_class_combination_simple + 
                                   arriv_diff_z + (1|position) + (1|dyad_ID) +
                                   (1|mm(initiator_JID, joiner_JID)), 
                                   data = dual_events,
                                   family = Gamma(link = "log"),
                                   prior = c(
                                     prior(normal(log(35), 0.5), class = "Intercept"),
                                     prior(normal(0, 0.3), class = "b"),
                                     prior(exponential(1), class = "sd"),
                                     prior(gamma(2, 0.5), class = "shape")
                                   ),
                                   iter = 6000, warmup = 2000,
                                   save_pars = save_pars(all = TRUE)
)


#Subset data where initiator at primary perch, no censored response
#Full dataset, censored response
social_tolerance_brm2_gamma <- brm(primary_duration ~ relationship + arriv_diff_z + 
                                     age_class_combination_simple + (1|position) + (1|dyad_ID) +
                                     (1|mm(initiator_JID, joiner_JID)), 
                                   data = dual_events2,
                                   family = Gamma(link = "log"),
                                   prior = c(
                                     prior(normal(log(35), 0.5), class = "Intercept"),
                                     prior(normal(0, 0.3), class = "b"),
                                     prior(exponential(1), class = "sd"),
                                     prior(gamma(2, 0.5), class = "shape")
                                   ),
                                   save_pars = save_pars(all = TRUE)
)
#iter = 6000, warmup = 2000,
#control = list(adapt_delta = 0.99, max_treedepth = 15),


social_tolerance_brm2_lognormal <- brm(duration2 | cens(duration_cens) ~ relationship + arriv_diff +
                                         age_class_combination_simple + (1|position) + (1|dyad_ID) +
                                         (1|mm(initiator_JID, joiner_JID)),
                                       data = dual_events,
                                       family = lognormal(link = "identity"),
                                       prior = c(
                                         prior(normal(log(6), 0.5), class = "Intercept"),
                                         prior(normal(0, 0.3), class = "b"),
                                         prior(exponential(1), class = "sd"),
                                         prior(exponential(1), class = "sigma")
                                       ),
                                       iter = 4000, warmup = 2000,
                                       control = list(adapt_delta = 0.95, max_treedepth = 12),
                                       save_pars = save_pars(all = TRUE)
)

social_tolerance_brm2 <- social_tolerance_brm2_gamma

summary(social_tolerance_brm2)
summary(social_tolerance_brm2_gamma)
summary(social_tolerance_brm2_lognormal)

check_collinearity(social_tolerance_brm2)

social_tolerance_brm2_gamma_loo <- loo(social_tolerance_brm2_gamma, moment_match = TRUE)
social_tolerance_brm2_lognormal_loo <- loo(social_tolerance_brm2_lognormal, moment_match = TRUE)
loo_compare(social_tolerance_brm2_gamma_loo, social_tolerance_brm2_lognormal_loo)

#Posterior predictive checks
pp_check(social_tolerance_brm2_gamma, type = "dens_overlay", ndraws = 100) +
  scale_x_log10()

pp_check(social_tolerance_brm2, type = "dens_overlay", ndraws = 100) +
  xlim(0, 100)

pp_check(social_tolerance_brm2_gamma, type = "dens_overlay", ndraws = 100) +
  xlim(0, 100)

pp_check(social_tolerance_brm2_lognormal, type = "dens_overlay", ndraws = 100) +
  xlim(0, 100)

pp_check(social_tolerance_brm2, type = "stat", stat = "var", ndraws = 1000)
pp_check(social_tolerance_brm2, type = "bars", ndraws = 1000)
pp_check(social_tolerance_brm2, type = "hist", ndraws = 1000)

y <- model.frame(social_tolerance_brm2)[[1]]
yrep <- posterior_predict(social_tolerance_brm2)
# Mean variance ratio across posterior draws:
dispersion <- apply(yrep, 1, function(x) var(x - y))
mean(dispersion)

#Plot model
plot(social_tolerance_brm2)
launch_shinystan(social_tolerance_brm2)

#Posterior distribution
as_draws_df(social_tolerance_brm2)
mcmc_areas(social_tolerance_brm2)
mcmc_intervals(social_tolerance_brm2)

#Evaluation and interpretation
loo(social_tolerance_brm2, moment_match = TRUE)
fitted(social_tolerance_brm2, scale = "response")
conditional_effects(social_tolerance_brm2)
conditional_effects(social_tolerance_brm2,effects = "relationship")
bayes_R2(social_tolerance_brm2)
emmeans(social_tolerance_brm2, ~ relationship, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(social_tolerance_brm2)

#Extract predictions
fitted(social_tolerance_brm2)

hypothesis(social_tolerance_brm2, "relationshippair > 0")
hypothesis(social_tolerance_brm2, "relationshippair - relationshipkin > 0")
hypothesis(social_tolerance_brm2, "relationshippair - relationshipneighbour > 0")
hypothesis(social_tolerance_brm2, "relationshippair - relationshipsite resident > 0")
hypothesis(social_tolerance_brm2, "relationshippair - relationshipother site > 0")

hypothesis(social_tolerance_brm2, "relationshipkin > 0")
hypothesis(social_tolerance_brm2, "relationshipkin - relationshipneighbour > 0")
hypothesis(social_tolerance_brm2, "relationshipkin - relationshipsite resident > 0")
hypothesis(social_tolerance_brm2, "relationshipkin - relationshipother site > 0")

hypothesis(social_tolerance_brm2, "relationshipneighbour > 0")
hypothesis(social_tolerance_brm2, "relationshipneighbour - relationshipsite resident > 0")
hypothesis(social_tolerance_brm2, "relationshipneighbour - relationshipother site > 0")

hypothesis(social_tolerance_brm2, "relationshipsite resident < 0")
hypothesis(social_tolerance_brm2, "relationshipother site < 0")

hypothesis(social_tolerance_brm2, "relationshipsite resident - relationshipother site > 0")

#Posterior marginal 
dual_events$relationship <- factor(dual_events$relationship, levels = c("other site" , "site resident", "foraging associate", "neighbour", "kin", "pair"))

posterior_marginal <- dual_events %>%
  add_epred_draws(
    object = social_tolerance_brm2,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_relationship <- posterior_marginal %>%
  group_by(.draw, relationship) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_summary <- posterior_draw_relationship %>%
  group_by(relationship) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_relationship %>%
  compare_levels(epred, by = relationship)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = relationship,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary


#3.3. Forager benefit across 3 contexts ----

#Prepare data for specific analysis

table(visit_data_primary2$relationship2)
table(visit_data_primary2$pre_relationship2)

length(unique(visit_data_primary2$JID[visit_data_primary2$relationship2 == "bonded"]))

tapply(visit_data_primary3$visit_duration, visit_data_primary3$relationship2, mean)
tapply(visit_data_primary3$arriv_diff, visit_data_primary3$relationship2, mean)


#Bayesian model 

#Model just with prior
default_prior()

#Censoring durations detected as 0 

# Make a new column for modeling
visit_data_primary3$visit_duration <- as.numeric(visit_data_primary3$visit_duration)
visit_data_primary3$visit_duration2 <- visit_data_primary3$visit_duration

# Set censored values to an approximate detection limit
visit_data_primary3$visit_duration2[visit_data_primary3$visit_duration == 0] <- 1
summary(visit_data_primary3$visit_duration2)

# Create a censoring indicator
visit_data_primary3$visit_duration_cens <- ifelse(visit_data_primary3$visit_duration == 0, "left", "none")
table(visit_data_primary3$visit_duration_cens)

#Scale latency of arrival
visit_data_primary3$arriv_diff <- ifelse(is.na(visit_data_primary3$arriv_diff), visit_data_primary3$visit_duration, visit_data_primary3$arriv_diff)
visit_data_primary3$arriv_diff_z <- scale(visit_data_primary3$arriv_diff)

#Dyadic event duration (0 for solo events)
visit_data_primary3$duration <- ifelse(is.na(visit_data_primary3$duration), 0, visit_data_primary3$duration)
visit_data_primary3$duration_z <- scale(visit_data_primary3$duration)

#Subset data: only include events where the initiator is at the primary perch
dual_events2 <- subset(dual_events, initiator_perch == "primary")

#Full dataset, censored response
social_tolerance_brm3_prior <- brm(visit_duration2 | cens(visit_duration_cens) ~ relationship2 + arriv_diff_z + (1|position) +
                                   (1|mm(JID, next_JID)), 
                                   data = visit_data_primary3, 
                                   family = Gamma(link = "log"),
                                   prior = c(
                                     prior(normal(log(18), 0.5), class = "Intercept"),
                                     prior(normal(0, 0.3), class = "b"),
                                     prior(exponential(1), class = "sd"),
                                     prior(gamma(2, 0.5), class = "shape")),
                                   sample_prior = "only")

#Subset of data where initiator at primary perch, no censored response

#Prior predictive checks 
pp_check(social_tolerance_brm3_prior, ndraws = 100) +
  xlim(0, 100)

pp_check(social_tolerance_brm3_prior, type="dens_overlay", ndraws = 100) + 
  scale_x_log10()

yrep <- posterior_predict(social_tolerance_brm3_prior)
quantile(yrep, 0.85)  # 99.9th percentile
max(yrep)              # largest simulated count

#Model

#Full dataset, censored response
social_tolerance_brm3 <- brm(visit_duration2 | cens(visit_duration_cens) ~ relationship2 + (1|position) +
                                     (1|mm(JID, next_JID)), 
                                   data = visit_data_primary3, 
                                   family = Gamma(link = "log"),
                                   prior = c(
                                     prior(normal(log(18), 0.5), class = "Intercept"),
                                     prior(normal(0, 0.3), class = "b"),
                                     prior(exponential(1), class = "sd"),
                                     prior(gamma(2, 0.5), class = "shape"))
)

#Subset data where initiator at primary perch, no censored response
#Full dataset, censored response
#iter = 6000, warmup = 2000,
#control = list(adapt_delta = 0.99, max_treedepth = 15),

social_tolerance_brm3 <- social_tolerance_brm3_gamma

summary(social_tolerance_brm3)
summary(social_tolerance_brm3_gamma)
summary(social_tolerance_brm3_lognormal)

check_collinearity(social_tolerance_brm3)

social_tolerance_brm3_gamma_loo <- loo(social_tolerance_brm3_gamma, moment_match = TRUE)
social_tolerance_brm3_lognormal_loo <- loo(social_tolerance_brm3_lognormal, moment_match = TRUE)
loo_compare(social_tolerance_brm3_gamma_loo, social_tolerance_brm3_lognormal_loo)

#Posterior predictive checks
pp_check(social_tolerance_brm3, type = "dens_overlay", ndraws = 100) +
  scale_x_log10()

pp_check(social_tolerance_brm3, type = "dens_overlay", ndraws = 100) 

pp_check(social_tolerance_brm3, type = "stat", stat = "var", ndraws = 1000)
pp_check(social_tolerance_brm3, type = "bars", ndraws = 1000)
pp_check(social_tolerance_brm3, type = "hist", ndraws = 1000)

y <- model.frame(social_tolerance_brm3)[[1]]
yrep <- posterior_predict(social_tolerance_brm3)
# Mean variance ratio across posterior draws:
dispersion <- apply(yrep, 1, function(x) var(x - y))
mean(dispersion)

#Plot model
plot(social_tolerance_brm3)
launch_shinystan(social_tolerance_brm3)

#Posterior distribution
as_draws_df(social_tolerance_brm3)
mcmc_areas(social_tolerance_brm3)
mcmc_intervals(social_tolerance_brm3)

#Evaluation and interpretation
loo(social_tolerance_brm3, moment_match = TRUE)
fitted(social_tolerance_brm3, scale = "response")
conditional_effects(social_tolerance_brm3)
conditional_effects(social_tolerance_brm3,effects = "relationship2")
bayes_R2(social_tolerance_brm3)
emmeans(social_tolerance_brm3, ~ relationship2, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(social_tolerance_brm3)

#Extract predictions
fitted(social_tolerance_brm3)

hypothesis(social_tolerance_brm3, "relationship2other > 0")

hypothesis(social_tolerance_brm3, "relationshippair > 0")
hypothesis(social_tolerance_brm3, "relationshippair - relationshipkin > 0")
hypothesis(social_tolerance_brm3, "relationshippair - relationshipneighbour > 0")
hypothesis(social_tolerance_brm3, "relationshippair - relationshipsite resident > 0")
hypothesis(social_tolerance_brm3, "relationshippair - relationshipother site > 0")

hypothesis(social_tolerance_brm3, "relationshipkin > 0")
hypothesis(social_tolerance_brm3, "relationshipkin - relationshipneighbour > 0")
hypothesis(social_tolerance_brm3, "relationshipkin - relationshipsite resident > 0")
hypothesis(social_tolerance_brm3, "relationshipkin - relationshipother site > 0")

hypothesis(social_tolerance_brm3, "relationshipneighbour > 0")
hypothesis(social_tolerance_brm3, "relationshipneighbour - relationshipsite resident > 0")
hypothesis(social_tolerance_brm3, "relationshipneighbour - relationshipother site > 0")

hypothesis(social_tolerance_brm3, "relationshipsite resident < 0")
hypothesis(social_tolerance_brm3, "relationshipother site < 0")

hypothesis(social_tolerance_brm3, "relationshipsite resident - relationshipother site > 0")

social_tolerance_m3 <- glmmTMB(visit_duration + 0.001 ~ relationship2 + arriv_diff_z + (1|JID) + (1|next_JID), family = Gamma, data = visit_data_primary3)
summary(social_tolerance_m3)
Anova(social_tolerance_m3)


#3.4. Displacement ----

table(visit_data_primary4$displaced, visit_data_primary4$relationship2)
tapply(visit_data_primary4$time_diff_next, visit_data_primary4$relationship2, mean, na.rm = T)

#Model just with prior
default_prior()

social_tolerance_brm4_prior <- brm(displaced ~ relationship2 + (1|mm(JID, next_JID)) + 
                                   (1|position) + (1|dyad_ID), 
                                   data = visit_data_primary4, 
                                   family = bernoulli(link = "logit"),
                                   prior = c(
                                     prior(normal(0, 0.5), class = "b"),        
                                     prior(normal(-1.5, 0.5), class = "Intercept")),
                                   sample_prior = "only"
)

#Prior predictive checks 
pp_check(social_tolerance_brm4_prior, ndraws = 100)

#Model
social_tolerance_brm4 <- brm(displaced ~ relationship2 + age_class_combination + 
                                   (1|mm(JID, next_JID)) + (1|position) + (1|dyad_ID), 
                                   data = visit_data_primary4, 
                                   family = bernoulli(link = "logit"),
                                   prior = c(
                                     prior(normal(0, 0.5), class = "b"),        
                                     prior(normal(-1.5, 0.5), class = "Intercept")),
                                   control = list(adapt_delta = 0.85, max_treedepth = 10)
                             
)

#Summary
summary(social_tolerance_brm4)

#Check collinearity
check_collinearity(social_tolerance_brm4)

#Posterior predictive checks
pp_check(social_tolerance_brm4, type = "dens_overlay", ndraws = 100)

#Plot model
plot(social_tolerance_brm4)
launch_shinystan(social_tolerance_brm4)

#Posterior distribution
as_draws_df(social_tolerance_brm4)
mcmc_areas(social_tolerance_brm4)
mcmc_intervals(social_tolerance_brm4)

#Evaluation and interpretation
loo(social_tolerance_brm4)
fitted(social_tolerance_brm4, scale = "response")
conditional_effects(social_tolerance_brm4)
bayes_R2(social_tolerance_brm4)
emmeans(social_tolerance_brm4, ~ relationship2, type = "response")
emmeans(social_tolerance_brm4, ~ relationship2, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(social_tolerance_brm4)

#Extract predictions
fitted(social_tolerance_brm4)

#Sensitivity analysis and prior checks 
prior_summary(social_tolerance_brm4)

#Extract predictions
fitted(social_tolerance_brm4)

#Posterior marginal 
visit_data_primary4$relationship2 <- factor(visit_data_primary4$relationship2, levels = c("solo", "other", "bonded"))

#Reduced number of draws due to larger dataset (computationally intense)
posterior_marginal <- visit_data_primary4 %>%
  add_epred_draws(
    object = social_tolerance_brm4,
    re_formula = NA,
    ndraws = 1000
  )

#Posterior summary for text reporting 
posterior_draw_relationship <- posterior_marginal %>%
  group_by(.draw, relationship2) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_summary <- posterior_draw_relationship %>%
  group_by(relationship2) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_relationship %>%
  compare_levels(epred, by = relationship2)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = relationship2,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

hypothesis(social_tolerance_brm4, "relationship2other > 0")
hypothesis(social_tolerance_brm4, "relationship2solo > 0")
hypothesis(social_tolerance_brm4, "relationship2other - relationship2solo > 0")

#Halfeye plot of model predictions
displacement_halfeye_plot <- ggplot(posterior_draw_relationship, aes(y = relationship2, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "grey") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Probability of displacement", y = "Relationship") +
  xlim(c(0, 1))
displacement_halfeye_plot

displacement_violin_plot <- ggplot(visit_data_primary4, aes(x = relationship2, y = displaced)) +
  geom_jitter(
    aes(x = relationship2, y = displaced), 
    data = visit_data_primary4, 
    color = "black",
    size = 2,
    alpha = 0.05,
    inherit.aes = FALSE,
    width = 0.05,
    height = 0.05) +
  geom_violin(alpha = 1, position = position_dodge(width = 0.8), width = 0.7, fill = "white") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(size = 14, family = "Garamond")) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.8), color = "white") +
  labs(x = "Relationship", y = "Displacement") +
  scale_x_discrete(guide = guide_axis(angle = 15))
displacement_violin_plot


#3.5. Queueing individual ----

#Q: Can the queueing individual feed after the feeding bird has left

visit_data_primary5 <- subset(visit_data_primary2, !JID_pre_queuer_JID_same == "solo")
visit_data_primary5$JID_pre_queuer_JID_same <- ifelse(visit_data_primary5$JID_pre_queuer_JID_same == "same", 1, 0)
visit_data_primary5$pre_duration_z <- scale(visit_data_primary5$pre_duration)
visit_data_primary5$pre_arriv_diff_z <- scale(visit_data_primary5$pre_arriv_diff)

table(visit_data_primary5$JID_pre_queuer_JID_same)
2577 / (2935 + 2577)


#Model just with prior
default_prior()

social_tolerance_brm5_prior <- brm(JID_pre_queuer_JID_same ~ pre_relationship +
                                     pre_age_class_combination_simple + pre_arriv_diff_z +
                                     pre_duration_z + (1|JID) + (1|mm(pre_JID, pre_queuer_JID)) + 
                                     (1|pre_dyad_ID) + (1|position), 
                                   data = visit_data_primary5, 
                                   family = bernoulli(link = "logit"),
                                   prior = c(
                                     prior(normal(0, 0.5), class = "b"),        
                                     prior(normal(-0.1, 0.5), class = "Intercept")),
                                   sample_prior = "only"
)

#Prior predictive checks 
pp_check(social_tolerance_brm5_prior, ndraws = 100)

#Model
social_tolerance_brm5 <- brm(JID_pre_queuer_JID_same ~ pre_relationship + 
                             pre_age_class_combination_simple
                             + (1|JID) + (1|mm(pre_JID, pre_queuer_JID)) 
                             + (1|pre_dyad_ID) + (1|position), 
                             data = visit_data_primary5, 
                             family = bernoulli(link = "logit"),
                             prior = c(
                               prior(normal(0, 0.5), class = "b"),        
                               prior(normal(-0.1, 0.5), class = "Intercept")),
                             control = list(adapt_delta = 0.99, max_treedepth = 10),
                             iter = 6000, warmup = 2000
)

#pre_duration_z?

#Summary
summary(social_tolerance_brm5)

#Check collinearity
check_collinearity(social_tolerance_brm5)

#Posterior predictive checks
pp_check(social_tolerance_brm5, type = "dens_overlay", ndraws = 100)

#Plot model
plot(social_tolerance_brm5)
launch_shinystan(social_tolerance_brm5)

#Posterior distribution
as_draws_df(social_tolerance_brm5)
mcmc_areas(social_tolerance_brm5)
mcmc_intervals(social_tolerance_brm5)

#Evaluation and interpretation
loo(social_tolerance_brm5)
fitted(social_tolerance_brm5, scale = "response")
conditional_effects(social_tolerance_brm5)
conditional_effects(social_tolerance_brm5,effects = "pre_relationship")
bayes_R2(social_tolerance_brm5)
emmeans(social_tolerance_brm5, ~ pre_relationship, type = "response")
emmeans(social_tolerance_brm5, ~ pre_relationship, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(social_tolerance_brm5)

#Sensitivity analysis and prior checks 
prior_summary(social_tolerance_brm5)

#Extract predictions
fitted(social_tolerance_brm5)

#Posterior marginal 
visit_data_primary5$pre_relationship <- factor(visit_data_primary5$pre_relationship, levels = c("other site" , "site resident", "foraging associate", "neighbour", "kin", "pair"))

#Reduced number of draws due to larger dataset (computationally intense)
posterior_marginal <- visit_data_primary5 %>%
  add_epred_draws(
    object = social_tolerance_brm5,
    re_formula = NA
    )

#Posterior summary for text reporting 
posterior_draw_relationship <- posterior_marginal %>%
  group_by(.draw, pre_relationship) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_summary <- posterior_draw_relationship %>%
  group_by(pre_relationship) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_relationship %>%
  compare_levels(epred, by = pre_relationship)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = pre_relationship,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

#Pairwise contrasts for Supp Mat using emmeans
contrasts <- social_tolerance_brm5 %>%
  emmeans(~ relationship) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  median_qi()
contrasts

hypothesis(social_tolerance_brm5, "pre_relationshippair < 0")
hypothesis(social_tolerance_brm5, "pre_relationshippair - pre_relationshipkin < 0")
hypothesis(social_tolerance_brm5, "pre_relationshippair - pre_relationshipneighbour > 0")
hypothesis(social_tolerance_brm5, "pre_relationshippair - pre_relationshipsite resident > 0")
hypothesis(social_tolerance_brm5, "pre_relationshippair - pre_relationshipother site > 0")

hypothesis(social_tolerance_brm5, "pre_relationshipkin > 0")
hypothesis(social_tolerance_brm5, "pre_relationshipkin - pre_relationshipneighbour > 0")
hypothesis(social_tolerance_brm5, "pre_relationshipkin - pre_relationshipsite resident > 0")
hypothesis(social_tolerance_brm5, "pre_relationshipkin - pre_relationshipother site > 0")

hypothesis(social_tolerance_brm5, "pre_relationshipneighbour > 0")
hypothesis(social_tolerance_brm5, "pre_relationshipneighbour - pre_relationshipsite resident > 0")
hypothesis(social_tolerance_brm5, "pre_relationshipneighbour - pre_relationshipother site > 0")

hypothesis(social_tolerance_brm5, "pre_relationshipsite resident > 0")
hypothesis(social_tolerance_brm5, "pre_relationshipother site > 0")

hypothesis(social_tolerance_brm5, "pre_relationshipsite resident - pre_relationshipother site > 0")

#Halfeye plot of model predictions
queuer_halfeye_plot <- ggplot(posterior_draw_relationship, aes(y = pre_relationship, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "grey") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Probability that queueing individual next visitor", y = "Relationship") +
  xlim(c(0, 1))
queuer_halfeye_plot

social_tolerance_plot <- ggarrange(social_tolerance_violin_plot, social_tolerance_halfeye_plot, displacement_halfeye_plot, queuer_halfeye_plot, ncol = 2, nrow = 2,  widths = c(1, 1))
social_tolerance_plot

#4. Social learning ----  

#Scale numeric predictor
visit_data_primary_first_juv2$day_ringed_z <- scale(visit_data_primary_first_juv2$day_ringed)

#4.1. Day of first visit ----

#Model just with prior
default_prior()

options(contrasts = c("contr.sum", "contr.poly"))

social_learning_brm1_prior <- brm(day ~ dual_juv_parent + day_ringed_z +
                                  (1|JID) + (1|position),
                                  data = visit_data_primary_first_juv2, 
                                  family = Gamma(link = "log"),
                                  prior = c(
                                    prior(normal(log(200), 0.15), class = "Intercept"),  
                                    prior(normal(0, 0.5), class = "b"),
                                    prior(exponential(1), class = "sd"),
                                    prior(gamma(150, 1), class = "shape")),
                                  control = list(adapt_delta = 0.85, max_treedepth = 10),
                                  sample_prior = "only")


#Prior predictive checks 
pp_check(social_learning_brm1_prior, ndraws = 100)+
  xlim(50,350)

pp_check(social_learning_brm1_prior, ndraws = 100)+
  scale_x_log10()

yrep <- posterior_predict(social_learning_brm1_prior)
quantile(yrep, 0.97)  # 99.9th percentile
max(yrep)             # largest simulated count

#Model
social_learning_brm1 <- brm(day ~ dual_juv_parent + day_ringed_z +
                            (1|JID) + (1|position),
                            data = visit_data_primary_first_juv2, 
                            family = Gamma(link = "log"),
                            prior = c(
                              prior(normal(log(200), 0.15), class = "Intercept"),  
                              prior(normal(0, 0.5), class = "b"),
                              prior(exponential(1), class = "sd"),
                              prior(gamma(150, 1), class = "shape")),
)

#iter = 6000, warmup = 2000,
#control = list(adapt_delta = 0.99, max_treedepth = 15)

#Model summary
summary(social_learning_brm1)

#Check collinearity
check_collinearity(social_learning_brm1)

#Posterior predictive checks
pp_check(social_learning_brm1, type = "dens_overlay", ndraws = 100) 

pp_check(social_learning_brm1, type = "stat", stat = "var", ndraws = 1000)
pp_check(social_learning_brm1, type = "bars", ndraws = 1000)
pp_check(social_learning_brm1, type = "hist", ndraws = 1000)

y <- model.frame(social_learning_brm1)[[1]]
yrep <- posterior_predict(social_learning_brm1)
# Mean variance ratio across posterior draws:
dispersion <- apply(yrep, 1, function(x) var(x - y))
mean(dispersion)

#Plot model
plot(social_learning_brm1)
launch_shinystan(social_learning_brm1)

#Posterior distribution
as_draws_df(social_learning_brm1)
mcmc_areas(social_learning_brm1)
mcmc_intervals(social_learning_brm1)

#Evaluation and interpretation
loo(social_learning_brm1, moment_match = TRUE)
fitted(social_learning_brm1, scale = "response")
conditional_effects(social_learning_brm1)
bayes_R2(social_learning_brm1)
emmeans(social_learning_brm1, ~ dual_juv_parent, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(social_learning_brm1)

#Extract predictions
fitted(social_learning_brm1)

#Posterior marginal 
visit_data_primary_first_juv2$dual_juv_parent <- factor(visit_data_primary_first_juv2$dual_juv_parent, levels = c("2", "1", "0"))

posterior_marginal <- visit_data_primary_first_juv2 %>%
  add_epred_draws(
    object = social_learning_brm1,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_relationship <- posterior_marginal %>%
  group_by(.draw, dual_juv_parent) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_summary <- posterior_draw_relationship %>%
  group_by(dual_juv_parent) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_relationship %>%
  compare_levels(epred, by = dual_juv_parent)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = dual_juv_parent,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

#Pairwise contrasts for Supp Mat using emmeans
contrasts <- social_learning_brm1 %>%
  emmeans(~ dual_juv_parent) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  median_qi()
contrasts

hypothesis(social_learning_brm1, "dual_juv_parent1 < 0")
hypothesis(social_learning_brm1, "dual_juv_parent2 < 0")
hypothesis(social_learning_brm1, "dual_juv_parent2 - dual_juv_parent1 < 0")

#Social learning halfeye plot 1 (day of first visit)
#Halfeye plot of model predictions
social_learning_halfeye_plot1 <- ggplot(posterior_draw_relationship, aes(y = dual_juv_parent, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "grey") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Day of first visit", y = "Following status") + 
  scale_y_discrete(labels=c("followed", "sometimes independent", "always independent"))
social_learning_halfeye_plot1

#Social learning violin plot 1 (day of first visit)
visit_data_primary_first_juv2$dual_juv_parent <- factor(visit_data_primary_first_juv2$dual_juv_parent, levels = c("0","1","2"))

social_learning_violin_plot1 <- ggplot(visit_data_primary_first_juv2, aes(x = dual_juv_parent, y = day)) +
  geom_violin(alpha = 1, position = position_dodge(width = 0.8), width = 0.7, fill = "white") +
  geom_jitter(
    aes(x = dual_juv_parent, y = day), 
    data = visit_data_primary_first_juv2, 
    color = "black",
    size = 2,
    alpha = 0.1,
    inherit.aes = FALSE,
    width = 0.2) +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(size = 14, family = "Garamond")) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.8), color = "black") +
  labs(x = "Following status", y = "Day of first visit") +
  scale_x_discrete(labels=c("always independent", "sometimes independent", "followed"))
social_learning_violin_plot1

#Panel 
social_learning_plot1 <- ggarrange(social_learning_violin_plot1, social_learning_halfeye_plot1, ncol = 2, nrow = 1,  widths = c(1, 1))
social_learning_plot1


#4.2. Number of visits ----

#Model on visit number by juveniles
tapply(visit_data_primary_first_juv2$visit_no, visit_data_primary_first_juv2$dual_juv_parent, mean)
tapply(visit_data_primary_first_juv2$visit_no, visit_data_primary_first_juv2$dual_juv_parent, median)

visit_data_primary_first_juv2$dual_juv_parent <- factor(visit_data_primary_first_juv2$dual_juv_parent, levels = c("0", "1", "2"))

visit_data_primary_first_juv2$obs <- row.names(visit_data_primary_first_juv2)

#Model just with prior
default_prior()

options(contrasts = c("contr.sum", "contr.poly"))

social_learning_brm2_prior <- brm(visit_no | trunc(lb = 1) ~ dual_juv_parent + 
                                  day_ringed_z + (1|JID) + (1|position),
                                  data = visit_data_primary_first_juv2, 
                                  family = poisson(link = "log"),
                                  prior = c(
                                    prior(normal(log(9), 4), class = "Intercept"),  
                                    prior(normal(0, 0.5), class = "b"),
                                    prior(exponential(1), class = "sd")),
                                  control = list(adapt_delta = 0.90),
                                  sample_prior = "only")

#iter = 6000, warmup = 2000,
#control = list(adapt_delta = 0.99, max_treedepth = 15),

#Prior predictive checks 
pp_check(social_learning_brm2_prior, ndraws = 100)

pp_check(social_learning_brm2_prior, ndraws = 100)+
  scale_x_log10()

yrep <- posterior_predict(social_learning_brm2_prior)
quantile(yrep, 0.97)  # 99.9th percentile
max(yrep)             # largest simulated count

#Model
social_learning_brm2 <- brm(visit_no | trunc(lb = 1) ~ dual_juv_parent + 
                            day_ringed_z + (1|JID) + (1|position),
                            data = visit_data_primary_first_juv2, 
                            family = poisson(link = "log"),
                            prior = c(
                              prior(normal(log(9), 4), class = "Intercept"),  
                              prior(normal(0, 0.5), class = "b"),
                              prior(exponential(1), class = "sd")),
                            control = list(adapt_delta = 0.90, max_treedepth = 15),
                            iter = 6000, warmup = 2000
)

social_learning_brm2_obs <- brm(visit_no | trunc(lb = 1) ~ dual_juv_parent + 
                              day_ringed_z + (1|JID) + (1|position) + (1|obs),
                            data = visit_data_primary_first_juv2, 
                            family = poisson(link = "log"),
                            prior = c(
                              prior(normal(log(9), 4), class = "Intercept"),  
                              prior(normal(0, 0.5), class = "b"),
                              prior(exponential(1), class = "sd")),
                            control = list(adapt_delta = 0.90, max_treedepth = 15),
                            iter = 6000, warmup = 2000
)


#iter = 6000, warmup = 2000,
#control = list(adapt_delta = 0.99, max_treedepth = 15)

#Model summary
summary(social_learning_brm2)

#Check collinearity
check_collinearity(social_learning_brm2)

#Posterior predictive checks
pp_check(social_learning_brm2, type = "dens_overlay", ndraws = 100) 
  
pp_check(social_learning_brm2, type = "dens_overlay", ndraws = 100) +
  scale_x_log10()

pp_check(social_learning_brm2_obs, type = "dens_overlay", ndraws = 100) +
  scale_x_log10()

pp_check(social_learning_brm2, type = "stat", stat = "var", ndraws = 1000)
pp_check(social_learning_brm2, type = "bars", ndraws = 1000)
pp_check(social_learning_brm2, type = "hist", ndraws = 1000)

y <- model.frame(social_learning_brm2)[[1]]
yrep <- posterior_predict(social_learning_brm2)
# Mean variance ratio across posterior draws:
dispersion <- apply(yrep, 1, function(x) var(x - y))
mean(dispersion)

#Plot model
plot(social_learning_brm2)
launch_shinystan(social_learning_brm2)

#Posterior distribution
as_draws_df(social_learning_brm2)
mcmc_areas(social_learning_brm2)
mcmc_intervals(social_learning_brm2)

#Evaluation and interpretation
loo(social_learning_brm2, moment_match = TRUE)
fitted(social_learning_brm2, scale = "response")
conditional_effects(social_learning_brm2)
bayes_R2(social_learning_brm2)
emmeans(social_learning_brm2, ~ dual_juv_parent, type = "response") 
emmeans(social_learning_brm2, ~ dual_juv_parent, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(social_learning_brm2)

#Extract predictions
fitted(social_learning_brm2)

#Posterior marginal 
visit_data_primary_first_juv2$dual_juv_parent <- factor(visit_data_primary_first_juv2$dual_juv_parent, levels = c("2", "1", "0"))

posterior_marginal <- visit_data_primary_first_juv2 %>%
  add_epred_draws(
    object = social_learning_brm2,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_relationship <- posterior_marginal %>%
  group_by(.draw, dual_juv_parent) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_summary <- posterior_draw_relationship %>%
  group_by(dual_juv_parent) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_relationship %>%
  compare_levels(epred, by = dual_juv_parent)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = dual_juv_parent,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

#Pairwise contrasts for Supp Mat using emmeans
contrasts <- social_learning_brm2 %>%
  emmeans(~ dual_juv_parent) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  median_qi()
contrasts

hypothesis(social_learning_brm2, "dual_juv_parent1 < 0")
hypothesis(social_learning_brm2, "dual_juv_parent2 < 0")
hypothesis(social_learning_brm2, "dual_juv_parent2 - dual_juv_parent1 < 0")

#Social learning halfeye plot 2 (number of visits)
#Halfeye plot of model predictions
social_learning_halfeye_plot2 <- ggplot(posterior_draw_relationship, aes(y = dual_juv_parent, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "grey") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Visit number", y = "Following status") +
  xlim(c(0, 25)) +
  scale_y_discrete(labels=c("followed", "sometimes independent", "always independent"))
social_learning_halfeye_plot2

#Social learning violin plot 2 (visit number)
visit_data_primary_first_juv2$dual_juv_parent <- factor(visit_data_primary_first_juv2$dual_juv_parent, levels = c("0","1","2"))

social_learning_violin_plot2 <- ggplot(visit_data_primary_first_juv2, aes(x = dual_juv_parent, y = visit_no)) +
  geom_violin(alpha = 1, position = position_dodge(width = 0.8), width = 0.7, fill = "white") +
  geom_jitter(
    aes(x = dual_juv_parent, y = visit_no), 
    data = visit_data_primary_first_juv2, 
    color = "black",
    size = 2,
    alpha = 0.1,
    inherit.aes = FALSE,
    width = 0.2) +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(size = 14, family = "Garamond")) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.8), color = "black") +
  labs(x = "Following status", y = "Visit number") +
  scale_x_discrete(labels=c("always independent", "sometimes independent", "followed")) + 
  ylim(0, 200) 
social_learning_violin_plot2

#Panel 
social_learning_plot2 <- ggarrange(social_learning_violin_plot2, social_learning_halfeye_plot2, ncol = 2, nrow = 1,  widths = c(1, 1))
social_learning_plot2

#Social learning overall panel plot
social_learning_plot <- ggarrange(social_learning_violin_plot1, social_learning_halfeye_plot1, social_learning_violin_plot2, social_learning_halfeye_plot2, ncol = 2, nrow = 2,  widths = c(1, 1))
social_learning_plot

#How many of the "following events" were also "dyadic queueing events"? 
visit_data_primary_first_juv$pre_dual_event_key <- visit_data_primary$pre_dual_event_key[match(visit_data_primary_first_juv$event_id, visit_data_primary$event_id)]
table(visit_data_primary_first_juv$parent, visit_data_primary_first_juv$pre_dual_event_key)

visit_data_primary_first_juv2$pre_dual_event_key <- visit_data_primary$pre_dual_event_key[match(visit_data_primary_first_juv2$event_id, visit_data_primary$event_id)]
table(visit_data_primary_first_juv2$parent, visit_data_primary_first_juv2$pre_dual_event_key)


#(4) SOCIAL NETWORK ----

#Tolerance/proximity networks

#dual events initiator and joiner ID
dual_events$visit_number_dyad <- dual_events_dyads$visit_number_dyad[match(dual_events$dyad_ID, dual_events_dyads$dyad_ID)]

dual_events_network <- dual_events[,c("joiner_JID", "initiator_JID", "dyad_ID", "relationship", "visit_number_dyad")]

#dual_event_network_dyads <- as.data.frame(table(dual_events_network$initiator_joiner_ID))
#dual_event_network_dyads$initiator_joiner_ID <- dual_event_network_dyads$Var1
#dual_event_network_dyads$visit_number_ini_joi <- dual_event_network_dyads$Freq
#dual_event_network_dyads <- subset(dual_event_network_dyads, select = -c(Var1, Freq))

#dual_events_network$visit_number_ini_joi <- dual_event_network_dyads$visit_number_ini_joi[match(dual_events_network$initiator_joiner_ID, dual_event_network_dyads$initiator_joiner_ID)]
#dual_events_network$weight <-dual_events_network$visit_number_ini_joi

dual_events_network$weight <- dual_events_network$visit_number_dyad

#dual_events_network <- dual_events_network %>% arrange(initiator_joiner_ID)

dual_events_network <- dual_events_network %>% 
  distinct(dyad_ID, .keep_all = T)

#Graph from association data

#Not directed
g_dual <- graph_from_data_frame(dual_events_network, directed=F, vertices= dual_events_individuals)

#Directed
g_dual <- graph_from_data_frame(dual_events_network, directed=T, vertices= dual_events_individuals)

#get edge list
get.edgelist(g_dual)

#get edge attribute
get.edge.attribute(g_dual)

#network characteristics 
sum(am_dual >0) #number of edges 1460
mean(rowSums(am_dual)) #mean number of associates: 17.44
mean(am_dual[am_dual >0]) #mean edge weight: 2.49

#vertex attributes
g_dual <- set_vertex_attr(g_dual, "sex", value = dual_events_individuals$sex)
#g_dual <- set_vertex_attr(g_dual , "age", value = dual_events_individuals$age24)
#g_dual <- set_vertex_attr(g_dual, "pref_site", value = dual_events_individuals$pref_site)
#g_dual <- set_vertex_attr(g_dual , "pair_ID", value = dual_events_individuals$pair_ID)

get.vertex.attribute(g_dual)

#subset graph
edge_weight_dual <- E(g_dual)$weight
edge_weight_duals <- E(g_duals)$weight

hist(E(g_dual)$weight, breaks = 100)

quantile(edge_weight_dual, probs = seq(0, 1, 1/100))
quantile(edge_weight_duals, probs = seq(0, 1, 1/100))

edges_to_keep_ids <- which(E(g_dual)$weight > 3)
g_duals <- subgraph.edges(graph = g_dual, eids = edges_to_keep_ids, delete.vertices = FALSE)

g_duals2 <- delete_vertices(g_duals, degree(g_duals)==0)

#colours
E(g_dual)$colour <-  E(g_dual)$relationship
E(g_duals)$colour <- E(g_duals)$relationship
E(g_duals2)$colour <- E(g_duals2)$relationship

#V(g_dual)$colour <- ifelse(V(g_dual)$pref_site == "Y", "#009E73", "#E69F00")
V(g_dual)$colour <- ifelse(V(g_dual)$sex == "M",  "#009E73", ifelse(V(g_dual)$sex == "F", "#E69F00", "#FFFFFF"))
#V(g_duals)$colour <- ifelse(V(g_duals)$pref_site == "Y", "#009E73", "#E69F00")
V(g_duals)$colour <- ifelse(V(g_duals)$sex == "M",  "#009E73", ifelse(V(g_dual)$sex == "F", "#E69F00", "#FFFFFF"))
#V(g_duals2)$colour <- ifelse(V(g_duals2)$pref_site == "Y", "#009E73", "#E69F00")
V(g_duals2)$colour <- ifelse(V(g_duals2)$sex == "M",  "#009E73", ifelse(V(g_dual)$sex == "F", "#E69F00", "#FFFFFF"))

E(g_dual)$colour <- ifelse(E(g_dual)$relationship == "pair", "#661100", ifelse(E(g_dual)$relationship == "kin", "#332288", "grey"))
E(g_duals)$colour <- ifelse(E(g_duals)$relationship == "pair", "#661100", ifelse(E(g_duals)$relationship == "kin", "#332288", "grey"))
E(g_duals)$colour <- ifelse(E(g_duals)$relationship == "pair", "#661100", "grey")
E(g_duals2)$colour <- ifelse(E(g_duals2)$relationship == "pair", "#661100", ifelse(E(g_duals2)$relationship == "kin", "#332288", "grey"))
E(g_duals2)$colour <- ifelse(E(g_duals2)$relationship == "pair", "#661100", "grey")

female_nodes <- which(V(g_dual)$sex=="F") 
male_nodes <- which(V(g_dual)$sex=="M") 
nonsexed_nodes <- which(V(g_dual)$sex== "U")

V(g_dual)$shape[female_nodes] <- "circle"
V(g_dual)$shape[male_nodes] <- "circle"
V(g_dual)$shape[nonsexed_nodes] <- "circle"

gorder(g_duals)

coords <- layout_(g_duals, as_star())
coords <- layout_(g_duals, in_circle())
coords <- layout_(g_duals, as_tree())
coords <- layout_(g_duals, nicely())
coords <- layout_(g_duals, on_grid())
coords <- layout_(g_duals, on_sphere())
coords <- layout_(g_duals, randomly())
coords <- layout_(g_duals, with_dh())
coords <- layout_(g_duals, with_fr())
coords <- layout_(g_duals, with_gem())
coords <- layout_(g_duals, with_graphopt())
coords <- layout_(g_duals, with_kk())
coords <- layout_(g_duals, with_lgl())
coords <- layout_(g_duals, with_mds())
coords <- layout_(g_duals, with_sugiyama())
coords <- layout_(g_duals, merge_coords())
coords <- layout_(g_duals, norm_coords())
coords <- layout_(g_duals, normalize())

coords <- layout_(g_duals, on_sphere())
coords <- layout_(g_duals2, on_sphere())

plot(g_dual, layout = coords, vertex.size = 4, vertex.label = NA, vertex.color = V(g_dual)$colour, edge.color = E(g_dual)$colour, edge.width = E(g_dual)$weight * 0.2, arrow.width = 0.25, edge.arrow.size= 0.25, edge.curved = 0.25)

par(mfrow=c(1,1))

windowsFonts(A = windowsFont("Garamond"))  

plot(g_duals, layout = coords, vertex.size = 4, vertex.label = NA, vertex.color = V(g_duals)$colour, edge.color = E(g_duals)$colour, edge.width = E(g_duals)$weight * 0.25, arrow.width = 0.25, edge.arrow.size= 0.25,  edge.curved = 0.25, family = "A")
plot(g_duals, layout = coords, vertex.size = 4, vertex.label = NA, vertex.color = V(g_duals)$colour, edge.color = E(g_duals)$colour, edge.width = E(g_duals)$weight * 0.25, arrow.width = 0.25, edge.arrow.size= 0.35,  edge.curved = 0.25, family = "A")
plot(g_duals, layout = coords, vertex.size = 4, vertex.label = NA, vertex.color = V(g_duals)$colour, edge.color = E(g_duals)$colour, edge.width = E(g_duals)$weight * 0.15,  edge.curved = 0.25, family = "A")

plot(g_duals2, layout = coords, vertex.size = 4, vertex.label = NA, vertex.color = V(g_duals2)$colour, edge.color = E(g_duals2)$colour, edge.width = E(g_duals2)$weight * 0.25, arrow.width = 0.25, edge.arrow.size= 0.25,  edge.curved = 0.25)

legend(x= 1.2, y= -0.1, c("Female","Male", "Unsexed"), pch=21, pt.bg = c("#E69F00", "#009E73", "white"), pt.cex= 1.2, cex= 1.2, bty="n", ncol=1)
legend(x= 1.2, y= -0.5, c("Kin","Pair", "Other"), lty = 1,  lwd = 3, col = c("#332288", "#661100", "grey"), cex= 1.2, bty="n", ncol=1)

#Centrality 
degree <- degree(g)
indegree <- degree(g, mode="in")
outdegree <- degree(g, mode="out")

degree <- as.data.frame(cbind(degree, indegree, outdegree))
degree <- cbind(JID = rownames(degree), degree)
rownames(degree) <- 1:nrow(degree)

degree_visits <- merge(x = degree, y = visit_duration, by = "JID", all = TRUE)
degree_visits <- merge(x = degree_visits, y = visits, by = "JID", all = TRUE)

degree_visits <- merge(x = degree, y = dual_event_duration, by = "JID", all = TRUE)
degree_visits <- merge(x = degree_visits, y = dual_event_visits, by = "JID", all = TRUE)

degree_visits$degree_per_visit <- (degree_visits$degree / degree_visits$visit_number)

eigen <- as.data.frame(eigen_centrality(g, directed = TRUE, scale = TRUE, weights = NULL))
eigen <- as.data.frame(eigen_centrality(g, directed = FALSE, scale = TRUE, weights = NULL))

eigen <- cbind(JID = rownames(eigen), eigen)
rownames(eigen) <- 1:nrow(eigen)
eigen <- subset(eigen[,1:2])

between <- as.data.frame(betweenness(g, directed = TRUE, weights = NULL))
between <- as.data.frame(betweenness(g, directed = FALSE, weights = NULL))
between <- betweenness(g, directed = FALSE, weights = NULL)

between <- cbind(JID = rownames(between), between)
rownames(between) <- 1:nrow(between)
between$between <- between[,2]
between <- subset(between[,c(1,3)])

centrality_visits <- merge(x = degree_visits, y = eigen, by = "JID", all = TRUE)
centrality_visits <- merge(x = centrality_visits , y = between, by = "JID", all = TRUE)

centrality_measures <- subset(centrality_visits[,c(2,8,9)])
centrality_measures <- na.omit(centrality_measures)
centrality_matrix = cor(centrality_measures)
centrality_matrix
corrplot(centrality_matrix)

plot(centrality_visits$degree, centrality_visits$vector)
plot(centrality_visits$degree, centrality_visits$visit_duration)
plot(centrality_visits$indegree, centrality_visits$visit_duration)
plot(centrality_visits$outdegree, centrality_visits$visit_duration)
plot(centrality_visits$degree_per_visit, centrality_visits$visit_duration)
plot(centrality_visits$vector, centrality_visits$visit_duration)

plot(centrality_visits$degree, centrality_visits$visit_number)
plot(centrality_visits$vector, centrality_visits$visit_number)
plot(centrality_visits$between, centrality_visits$visit_number)

cor.test(centrality_visits$degree, centrality_visits$visit_number)
cor.test(centrality_visits$vector, centrality_visits$visit_number)
cor.test(centrality_visits$between, centrality_visits$visit_number)
cor(c(centrality_visits$degree, centrality_visits$vector, centrality_visits$between))

boxplot(dual_events_individuals$degree ~ dual_events_individuals$sex)
plot(dual_events_individuals$degree ~ dual_events_individuals$age)

boxplot(dual_events_individuals$degree ~ dual_events_individuals$sex)
boxplot(dual_events_individuals$eigenvector ~ dual_events_individuals$sex)
boxplot(dual_events_individuals$between~ dual_events_individuals$sex)

plot(dual_events_individuals$degree ~ dual_events_individuals$age)
plot(dual_events_individuals$eigenvector ~ dual_events_individuals$age)
plot(dual_events_individuals$between~ dual_events_individuals$age)
plot(dual_events$initiator_degree ~ dual_events$joiner_degree)

#Distance between nodes based on edge weight
plot(g_duals2, layout=layout.fruchterman.reingold(g_duals2,weights=E(g_duals2)$weight^3), vertex.size = 2.5, vertex.label = NA, vertex.color = "black", arrow.width = 0.1, edge.arrow.size= 0.1, edge.curved = 0.25) 

community <- leading.eigenvector.community(g)
communities <- list()
for (i in 1:max(community$membership)) {
  communities[[i]] <- which(community$membership == i)
}

plot(g, vertex.color= V(g)$colour, edge.arrow.size= 0.05, arrow.width = 0.05, edge.width=E(g)$weight, vertex.size= 4, vertex.label=NA)
plot(g, vertex.color= V(g)$colour, vertex.color = "circle", edge.width=E(g)$weight, vertex.size= 4, vertex.label=NA)

par(mfrow=c(1,1))
par(mfrow=c(2,2), mar=c(1,1,1,1))
plot(g, layout=layout.sphere, main="sphere")
plot(g, layout=layout.circle, main="circle")
plot(g, layout=layout.random, main="random")
plot(g, layout=layout.fruchterman.reingold, main="fruchterman.reingold")

#(5) MAP ----
feedercoord <- read.csv("feedercoord.csv", header = T, stringsAsFactors = F)

feedercoordY <- subset(feedercoord, feedercoord$site == "Y")
feedercoordZ <- subset(feedercoord, feedercoord$site == "Z")

library(ggmap)
library(ggsn)
library(ggrepel)
register_google(key = "AIzaSyDIvL-XLO1etdpwxUOfviz0wJ84n_et3rA")

install.packages("extrafont")
library(extrafont)
font_import()
loadfonts()

#Map for both sites
mean(feedercoord$latitude)
mean(feedercoord$longitude)
location <- c(lon = -5.176987, lat = 50.19348)

pointlabels <-annotate("text", x = feedercoord$longitude,y= feedercoord$latitude,size=3,label = as.vector(feedercoord$position))

feedermap <- get_map(location=location,
                     source="google", maptype="satellite", crop=FALSE, zoom = 15, color = "bw")

ggmap(feedermap)

#Dyadic- and single-perch feeding stations
feedermap <- ggmap(feedermap) +
  geom_point(aes(x = longitude, y = latitude, fill = feederfun), data = feedercoord, 
             alpha = 1, size = 4, shape = 21, color = "white") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none", text = element_text(size = 20, family = "Garamond")) +
  scale_fill_manual(values = c("solo" = "darkgrey", "dual" = "white", "both" = "white")) +
  labs(color='Feeder function') +
  xlab("Longitude (W)") +
  ylab("Latitude (N)") +
  scalebar(location = "bottomright", x.min = -5.175, x.max = -5.165,
           y.min = 50.187, y.max = 50.192,
           dist = 250, transform = TRUE, dist_unit = "m", model = 'WGS84',
           box.fill = c("black", "white"), st.color = 'white', box.color = 'white', 
           height = 0.05, st.dist = 0.1, border.size = 1)  +
  geom_label_repel(
    aes(x = longitude, y = latitude, label = position),
    data= feedercoord,
    size = 4, 
    box.padding = 0.2, point.padding = 0.3,
    segment.color = 'grey50')
feedermap

#Dyadic feeding stations only
feedercoord_dyadic <- subset(feedercoord, !feedercoord$feederfun == "solo")

feedermap <- ggmap(feedermap) +
  geom_point(aes(x = longitude, y = latitude, fill = feederfun), data = feedercoord_dyadic, 
             alpha = 1, size = 4, shape = 21, color = "white") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none", text = element_text(size = 20, family = "Garamond")) +
  scale_fill_manual(values = c("solo" = "darkgrey", "dual" = "white", "both" = "white")) +
  labs(color='Feeder function') +
  xlab("Longitude (W)") +
  ylab("Latitude (N)") +
  scalebar(location = "bottomright", x.min = -5.175, x.max = -5.165,
           y.min = 50.187, y.max = 50.192,
           dist = 250, transform = TRUE, dist_unit = "m", model = 'WGS84',
           box.fill = c("black", "white"), st.color = 'white', box.color = 'white', 
           height = 0.05, st.dist = 0.1, border.size = 1)  +
  geom_label_repel(
    aes(x = longitude, y = latitude, label = position),
    data= feedercoord_dyadic,
    size = 4, 
    box.padding = 0.2, point.padding = 0.3,
    segment.color = 'grey50')
feedermap



#Archive 
feedermap <- ggmap(feedermap) +
  geom_point(aes(x = longitude, y = latitude, color = feederfun), data = feedercoord, 
             alpha = 1, size = 4) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("solo" = "yellow", "dual" = "#FF6633", "both" = "yellow")) +
  labs(color='Feeder function') +
  xlab("Longitude (W)") +
  ylab("Latitude (N)") +
  scalebar(location = "bottomright", x.min = -5.175, x.max = -5.165,
           y.min = 50.187, y.max = 50.192,
           dist = 250, transform = TRUE, dist_unit = "m", model = 'WGS84',
           box.fill = c("black", "white"), st.color = 'white', box.color = 'white', 
           height = 0.05, st.dist = 0.1, border.size = 1) 
feedermap

feedermap <- ggmap(feedermap) +
  geom_point(aes(x = longitude, y = latitude, color = feederfun), data = feedercoord, 
             alpha = 1, size = 6) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("solo" = "white", "dual" = "white", "both" = "white")) +
  labs(color='Feeder function') +
  xlab("Longitude (W)") +
  ylab("Latitude (N)") +
  scalebar(location = "bottomright", x.min = -5.175, x.max = -5.165,
           y.min = 50.187, y.max = 50.192,
           dist = 250, transform = TRUE, dist_unit = "m", model = 'WGS84',
           box.fill = c("black", "white"), st.color = 'white', box.color = 'white', 
           height = 0.05, st.dist = 0.1, border.size = 1)  +
  pointlabels
feedermap

feedermap <- ggmap(feedermap) +
  geom_point(aes(x = longitude, y = latitude, fill = feederfun), data = feedercoord, 
             alpha = 1, size = 6, shape = 21, color = "white") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("solo" = "black", "dual" = "white", "both" = "white")) +
  labs(color='Feeder function') +
  xlab("Longitude (W)") +
  ylab("Latitude (N)") +
  scalebar(location = "bottomright", x.min = -5.175, x.max = -5.165,
           y.min = 50.187, y.max = 50.192,
           dist = 250, transform = TRUE, dist_unit = "m", model = 'WGS84',
           box.fill = c("black", "white"), st.color = 'white', box.color = 'white', 
           height = 0.05, st.dist = 0.1, border.size = 1)  +
  pointlabels
feedermap

#Map for Y
mean(feedercoordY$latitude)
mean(feedercoordY$longitude)
locationY <- c(lon = -5.182384, lat = 50.1905)

feedermapY <- get_map(location=locationY,
                      source="google", maptype="satellite", crop=FALSE, zoom = 17)

ggmap(feedermapY)

feedermapY <- ggmap(feedermapY) +
  geom_point(aes(x = longitude, y = latitude, color = feederfun), data = feedercoordY, 
             alpha = 1, size = 3) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("solo" = "yellow", "dual" = "red", "both" = "yellow")) +
  labs(color='Feeder function') +
  xlab("Longitude (W)") +
  ylab("Latitude (N)") +
  scalebar(location = "bottomright", x.min = -5.1795, x.max = -5.181,
           y.min = 50.1887, y.max = 50.1897,
           dist = 50, transform = TRUE, dist_unit = "m", model = 'WGS84',
           box.fill = c("black", "white"), st.color = 'white', box.color = 'white', 
           height = 0.1, st.dist = 0.1, border.size = 1) 
feedermapY
