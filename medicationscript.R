###### Medication quantification script
# Arjen Cupido
# 09-06-2021

# Load libraries
library(plyr)
library(tidyverse)

# Set working directory
setwd("C:/Users/arjen/Documents/research/UCLA/METSIM/data/R")

# Load data (dataset and medication identification list)

m <- readxl::read_xlsx('DIRECT_medicati.xlsx')
uniquemeds <- read.csv('uniquemeds.csv')

# Format from wide to long format, based on visit:

m0 <- m %>% 
  select(., c('METSIM ID', starts_with("D0"))) %>% 
  add_column(time_point="D0") %>% 
  rename_at(vars(starts_with("D0_")), list(~str_remove(., "D0_"))) %>%
  unite("MOMID_TP",  c("METSIM ID","time_point"), remove = FALSE)

m18 <- m %>% 
  select(., c('METSIM ID', starts_with("D18"))) %>% 
  add_column(time_point="D18") %>% 
  rename_at(vars(starts_with("D18_")), list(~str_remove(., "D18_"))) %>%
  unite("MOMID_TP",  c("METSIM ID","time_point"), remove = FALSE)

m48 <- m %>% 
  select(., c('METSIM ID', starts_with("D48"))) %>% 
  add_column(time_point="D48") %>% 
  rename_at(vars(starts_with("D48_")), list(~str_remove(., "D48_"))) %>%
  unite("MOMID_TP",  c("METSIM ID","time_point"), remove = FALSE)

# Combine data
d <- bind_rows(m0, m18, m48)

# format free text, get rid of special characters
d$medicati <- tolower(d$medicati)

d$medicati2 <- gsub('[[:digit:]]+', '', d$medicati)
d$medicati2 <- gsub(' mg', '', d$medicati2)
d$medicati2 <- gsub('/', '', d$medicati2)
d$medicati2 <- gsub('.g', '', d$medicati2)
d$medicati2 <- gsub(' ml', '', d$medicati2)

head(d$medicati2)

### Start of script identifying medication in free text:

classes <- unique(uniquemeds$class_indication)[-2] # remove NA class
classes

# Code medication

for(i in 1:length(classes)){

    # Select drug class
  x <- subset(uniquemeds, class_indication == classes[i]) # subset all meds in that specific class
  x <- paste(x$searchterm, collapse= "|", sep ="|") # Collapse the searchterms
  
  # Search in free text, write TRUE/FALSE in drug class variable
  d[paste0(classes[i])] <- grepl(x, d$medicati2)
}

# Write data including reference database.
writexl::write_xlsx(list(data=d, referencelist =uniquemeds), path = 'DIRECT_medicationlist2.xlsx')
