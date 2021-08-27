### METSIM DIETARY DATA
###packages
library("scales")
library("grid")
library("biom")
library("RColorBrewer")
library ("permute")
library ("lattice")
library("vegan")
library("DESeq2")
#library("ancom.R")
library("gtable")
library("grid")
library("gridExtra")
library("reshape2")
library("ggbeeswarm")
library("corrplot")
library('magrittr')
library('plyr')
library("phyloseq")
library("tidyverse") # ggplot2, tibble, tidyr, readr, purrr, dplyr, stringr, forcats
setwd("C:/Users/arjen/Documents/research/UCLA/METSIM/data/R")

`%notin%` <- negate(`%in%`)



###################### VARIABLES LIST ###############################

var_per <- c(#"MOMID",
             "DNA",
             "CityCode",
             "Sex",
             "Birthy",
             "Finnish",
             "Res_date",
             "B_CityCode",
             "B_ProvinceCode",
             "B_Weight",
             "B_Height")


var_env <- c("smoke",
             "smokeq",
             "smokey",
             "workex",
             "hobbyex",
             "F_smoke",
             "F_smokeq",
             "F_smokey",
             "F_workex",
             "F_hobbyex",
             "F2_smoke",
             "F2_smokeq",
             "F2_smokey",
             "F2_workex",
             "F2_hobbyex")

var_diet <- c("Freq_veg", 
              "Freq_fruit", 
              "Freq_leanfish",
              "Freq_fattyfish", 
              "Freq_shellfish", 
              "Freq_alclt3", 
              "Freq_alclt6", 
              "Freq_alcge6", 
              "Freq_wine", 
              "Freq_liqueur", 
              "Freq_strongwine", 
              "Freq_blend",
              "Q_alclt3", 
              "Q_alclt6",
              "Q_alcge6",
              "Q_wine",
              "Q_liqueur",
              "Q_strongwine",             
              "Q_blend",
              "Milk",           
              "Milk_quality",
              "Milk_quantity",
              "Dairy_other",
              "Spread_sat",
              "Spread_marg",
              "Spread_no",
              "Cookfat_sat",
              "Cookfat_oils",
              "Cookfat_marg",
              "Cookfat_no",
              "Redmeat_freq",
              "Redmeat_g",
              "Redmeat_gwk", 
              "Cheese_freq", 
              "Cheese_g",
              "Cheese_gvko",
              "Cheese_other",
              "Cereal_24_serv_wholegrain",
              "Cereal_24_serv_wheat",
              "Cereal_24_serv_pastry")



metadata <- read.csv("METSIM_metadata.csv", header = TRUE)
mapping <- read.table("METSIM_mapping.txt", header = TRUE)
NMRdata1 <-readxl:: read_xlsx('METSIM-DIRECT-D0_D18_2019-08-05.xlsx')
NMRdata2 <- readxl:: read_xlsx('METSIM-DIRECT_NMR_18710-30-Jun-2019-Results_2019-08-02.xlsx')

medication <- readxl::read_xlsx('DIRECT_medicationlist2.xlsx')


###################### DIRECT D0, D18, D48 DATA ##########################
#### Per Time Point
#1.Select correct time point
#2.Select only columns for that time point and the personal variables
#3.Add timepoint column
#4.Rename columns by removing timepoint designation so they are consistent between time points
#5.Make the date of sampling date (ResDate) common by removing the extra information 
#6.Create unique identifier that is MOMID_Timepoint
#7.Remove "BarcodeInQuestionnaire" because it doesn't allow the dataframes to join since time D48 becomes an integer
  #BC <- select(mapping, ends_with("BarcodeInQuestionnaire"))

##### create subset of each timepoint
D0 <- metadata %>% 
  filter(D0_MOMID != "NA") %>% 
  select(., c(var_per, starts_with("D0"))) %>% 
  add_column(time_point="D0") %>% 
  rename_at(vars(starts_with("D0_")), list(~str_remove(., "D0_"))) %>%
  rename_at(vars(ends_with("ResDate")), list(~"ResDate")) %>%
  unite("MOMID_TP",  c("MOMID","time_point"), remove = FALSE) %>%
  select(-"BarcodeInQuestionnaire")
  
D18 <- metadata %>% 
  filter(D18_MOMID != "NA") %>% 
  select(., c(var_per,starts_with("D18"))) %>% 
  add_column(time_point="D18") %>% 
  rename_at(vars(starts_with("D18_")), list(~str_remove(., "D18_"))) %>%
  rename_at(vars(ends_with("ResDate")), list(~"ResDate")) %>%
  unite("MOMID_TP",  c("MOMID","time_point"), remove = FALSE) %>%
  select(-"BarcodeInQuestionnaire")

D48 <- metadata %>% 
  filter(D48_MOMID != "NA") %>% 
  select(., c(var_per,starts_with("D48"))) %>% 
  add_column(time_point="D48") %>% 
  rename_at(vars(starts_with("D48_")), list(~str_remove(., "D48_"))) %>%
  rename_at(vars(ends_with("ResDate")), list(~"ResDate")) %>% 
  unite("MOMID_TP",  c("MOMID","time_point"), remove = FALSE) %>%
  select(-"BarcodeInQuestionnaire")

# NMR data

NMR0 <- NMRdata1 %>% 
  select(., c('MOMID', starts_with("D0"))) %>% 
  add_column(time_point="D0") %>% 
  rename_at(vars(starts_with("D0_")), list(~str_remove(., "D0_"))) %>%
  unite("MOMID_TP",  c("MOMID","time_point"), remove = FALSE) %>% 
  subset(!is.na(LDL_C)) # Remove people without values


# Timepoint zero data in the second datafile
NMR0_2 <- NMRdata2 %>% 
  select(., c('MOMID', starts_with("D0"))) %>% 
  add_column(time_point="D0") %>% 
  rename_at(vars(starts_with("D0_")), list(~str_remove(., "D0_"))) %>%
  unite("MOMID_TP",  c("MOMID","time_point"), remove = FALSE) %>% 
  subset(!is.na(LDL_C)) # Remove people without values

NMR0_2$MOMID_TP %in% NMR0$MOMID_TP # All samples in the second datafile can also be found in first datafile. 
NMR0_duplicates <- subset(NMR0, MOMID %in% NMR0_2$MOMID)
plot(x=NMR0_2$LDL_C, y = NMR0_duplicates$LDL_C) # datapoints don't overlap. Decided to keep data from the first dataset
writexl::write_xlsx(list(NMR0, NMR0_2), path="NMR0_data_incl_duplicates.xlsx")  # For reference


NMR18 <- NMRdata1 %>% 
  select(., c('MOMID', starts_with("D18"))) %>% 
  add_column(time_point="D18") %>% 
  rename_at(vars(starts_with("D18_")), list(~str_remove(., "D18_"))) %>%
  unite("MOMID_TP",  c("MOMID","time_point"), remove = FALSE) %>% 
  subset(!is.na(LDL_C))
NMR18_2 <- NMRdata2 %>% 
  select(., c('MOMID', starts_with("D18"))) %>% 
  add_column(time_point="D18") %>% 
  rename_at(vars(starts_with("D18_")), list(~str_remove(., "D18_"))) %>%
  unite("MOMID_TP",  c("MOMID","time_point"), remove = FALSE) %>% 
  rename_at(vars(ends_with("pct")), list(~str_remove(., "ct"))) %>%
  subset(!is.na(LDL_C))

NMR18full <- bind_rows(NMR18, NMR18_2)  
writexl::write_xlsx(list(NMR18, NMR18_2, NMR18full), path="NMR18_data.xlsx")  
NMR18_2$MOMID_TP %in% NMR18$MOMID_TP # None in earlier dataset. 
table(D18$MOMID_TP %in% NMR18full$MOMID_TP)



NMR48 <- NMRdata2 %>% 
  select(., c('MOMID', starts_with("D48"))) %>% 
  add_column(time_point="D48") %>% 
  rename_at(vars(starts_with("D48_")), list(~str_remove(., "D48_"))) %>%
  unite("MOMID_TP",  c("MOMID","time_point"), remove = FALSE) %>% 
  rename_at(vars(ends_with("pct")), list(~str_remove(., "ct"))) %>%
  subset(!is.na(LDL_C))
writexl::write_xlsx(list(NMR48), path="NMR48.xlsx")  
table(D48$MOMID_TP %in% NMR48$MOMID_TP)

###### Add together time point by rows

#need to know what columns are different because the have to be the same to join by row
# setdiff(D0,D18)
# setdiff(D0,D48)
# setdiff(D18,D48) #have all the same columns
# 
# #these are columns in D0 that aren't in D18/D48
# drop_cols <- c('fsdK01', 'P_lac30', 'S_tottg120', 'S_tottg90', 'S_tottg45', 'S_tottg15', 'P_lac45', 'S_tottg60', 'P_pyr120', 'P_pyr0', 'S_tottg30', 'P_lac60', 'fsdKGK', 'P_lac15', 'P_lac0', 'valueK01', 'valueLacBas', 'P_lac90', 'valueK12', 'fsdK12', 'valueKGK', 'P_pyr30', 'P_lac120')
# 
# #these are columns in D18/D48 that need to be added to D0 so all timepoints can be combined
# add_cols <- c('P_glucagon120', 'UAlbV', 'UAlbTh1', 'UAlbTm0', 'S_totalc', 'UAlbTm1', 'P_glucagon90', 'P_glucagon45', 'UAlbTh0', 'P_glucagon15', 'P_glucagon0', 'P_glucagon30', 'P_glucagon60', 'S_ldlc', 'S_hdlc')
# 
# #make D0 have the same columms as D18/D48 and D18/D48 the same as D0.
# D0_same <- D0 %>% 
#   select(-drop_cols) 
# for(x in add_cols){
#   D0_same[[x]] = NA
# }
# 
# setdiff(D0_same,D18)
# setdiff(D0_same,D48)
# 
# D01848 <- bind_rows(D0_same, D18, D48)
# D01848 <- left_join(mapping, D01848,  by = "MOMID_TP")
# 
# write.table(D01848, "METSIM_D_metadata.txt", quote = FALSE, sep='\t', row.names=F, col.names = T)

D01848full <- bind_rows(D0, D18, D48)
D01848microbiome <- left_join(mapping, D01848full, by = "MOMID_TP")

# NMR data
setdiff(NMR0, NMR18)
setdiff(NMR18, NMR18_2)
setdiff(NMR0, NMR48)
setdiff(NMR18, NMR48)

# Renaming on https://www.nature.com/articles/s41467-019-12703-7#MOESM1
NMR48 <- NMR48 %>% rename(bOHBut = bOHbutyrate,LDL_D = LDL_size, PC = Phosphatidylc,FAw3 = Omega_3, UnSat = Unsaturation,
                          DHA_FA = DHA_p, FAw3_FA = Omega_3_p,Serum_TG =Total_TG,PUFA_FA = PUFA_p,
                          Lac = Lactate, 
                          Crea = Creatinine,
                          Glc = Glucose, 
                          TG_PG = TG_by_PG , 
                          FAw6 = Omega_6, 
                          SM = Sphingomyelins, 
                          MUFA_FA = MUFA_p, 
                          FAw6_FA = Omega_6_p,
                          Cit = Citrate,
                          ApoB_ApoA1 = ApoB_by_ApoA1, 
                          Serum_C = Total_C,
                          SFA_FA = SFA_p,
                          Ace = Acetate,
                          AcAce = Acetoacetate,
                          TotPG = Phosphoglyc,
                          VLDL_D = VLDL_size,
                          TotCho = Cholines,
                          FreeC = Free_C,
                          TotFA = Total_FA,
                          LA_FA = LA_p,
                          Alb = Albumin,
                          EstC = Esterified_C,
                          HDL_D = HDL_size,
                          Gp = GlycA)
NMR18_2 <- NMR18_2 %>% rename(bOHBut = bOHbutyrate,LDL_D = LDL_size, PC = Phosphatidylc,FAw3 = Omega_3, UnSat = Unsaturation,
                              DHA_FA = DHA_p, FAw3_FA = Omega_3_p,Serum_TG =Total_TG,PUFA_FA = PUFA_p,
                              Lac = Lactate, 
                              Crea = Creatinine,
                              Glc = Glucose, 
                              TG_PG = TG_by_PG , 
                              FAw6 = Omega_6, 
                              SM = Sphingomyelins, 
                              MUFA_FA = MUFA_p, 
                              FAw6_FA = Omega_6_p,
                              Cit = Citrate,
                              ApoB_ApoA1 = ApoB_by_ApoA1, 
                              Serum_C = Total_C,
                              SFA_FA = SFA_p,
                              Ace = Acetate,
                              AcAce = Acetoacetate,
                              TotPG = Phosphoglyc,
                              VLDL_D = VLDL_size,
                              TotCho = Cholines,
                              FreeC = Free_C,
                              TotFA = Total_FA,
                              LA_FA = LA_p,
                              Alb = Albumin,
                              EstC = Esterified_C,
                              HDL_D = HDL_size,
                              Gp = GlycA)


NMR48$`Sample id`<-NULL
NMR18_2$`Sample id` <- NULL
#NMR0 <- NMR0 %>% select(-c(`Pyr`, `Gly`, `Glol`))
NMR01848 <- bind_rows(NMR0, NMR18, NMR18_2, NMR48)
sapply(NMR01848, function(x) sum(is.na(x)))
sapply(D01848microbiome, function(x) sum(is.na(x)))

# Write all NMR data
write.table(NMR01848, "METSIM_NMR_metadata.txt", quote = FALSE, sep='\t', row.names=F, col.names = T)

DNMR <- left_join(D01848microbiome, NMR01848, by = 'MOMID_TP')

# Write all data
write.csv(DNMR, "METSIM_all_data.csv")

# Subset all data to the ones we have microbiome on. 
DNMRmicrobiome <- left_join(D01848microbiome, NMR01848, by = 'MOMID_TP') 
DNMRmicrobiome <- subset(DNMRmicrobiome, DNMRmicrobiome$DNA != "NA")

sapply(DNMRmicrobiome, function(x) sum(is.na(x)))

# Write microbiome subset
write.csv(DNMRmicrobiome, "METSIM_microbiomesubset_NMR_Meta_data.csv")

missingpercentages <- data.frame(totalNA = sapply(DNMRmicrobiome, function(x) sum(is.na(x)))/1717*100) # Quite some missings. We'll remove the ones that are irrelevant
throwaway <- rownames(subset(missingpercentages, totalNA>=50)) # Throw away averything that is missing in > 50% of subjects. 

# Throw away other variables not deemed interesting
throwaway2<- c( 'MOMID.x', 'Folder', 'MOMFollowUpTime', 'time_point.y', 'MOMID.y', 'DNA', 'CityCode', 'Sex')

DNMRmicrobiomeclean <- DNMRmicrobiome %>% 
  select(-throwaway) %>% select(-throwaway2)
colnames(DNMRmicrobiomeclean[grep(colnames(DNMRmicrobiomeclean), pattern= ".y")])

View(data.frame(totalNA = sapply(DNMRmicrobiomeclean, function(x) sum(is.na(x)))/1717*100))


# Write cleaner dataset
write.csv(DNMRmicrobiomeclean, "CLEAN_METSIM_microbiomesubset_NMR_Meta_data.csv")

# Add medication data
d <- read.csv( "CLEAN_METSIM_microbiomesubset_NMR_Meta_data.csv")
m <- readxl::read_xlsx('DIRECT_medicationlist2.xlsx')[,-c(2,4,5)]
d <- left_join(d, m, by = "MOMID_TP")

write.csv(d, "CLEAN_METSIM_microbiomesubset_NMR_Meta_drug_data.csv")
d <- read.csv("CLEAN_METSIM_microbiomesubset_NMR_Meta_drug_data.csv")

# Add outcome and genetics 
metadata2 <- subset(metadata, select = c("MOMID",
                                         "TotalCADEvents_2013",
                                         "totalcw",
                                         "TotalT2DM_022017",
                                         "TotalT2DM_062016",
                                         "TotalT2DM_072017",
                                         "TotalT2DM_072018",
                                         "TotalT2DM_092014",
                                         "TotalT2DM_102018",
                                         "TotalT2DM_2013",
                                         "Reimb_CHD_date",
                                         "Reimb_CHD_date_122017",
                                         "Reimb_DM_date",
                                         "Reimb_DM_date_122017",
                                         "Reimb_HF_date",
                                         "Reimb_HF_date_122017",
                                         "Reimb_HT_date",
                                         "Reimb_HT_date_122017",
                                         "GRS_2hPG_9nw",
                                         "GRS_2hPG_9w",
                                         "GRS_BMI_95nw",
                                         "GRS_BMI_95w",
                                         "GRS_FPG_35nw",
                                         "GRS_FPG_35w",
                                         "GRS_HbA1c_11w",
                                         "GRS_HDL_67w",
                                         "GRS_IR_9nw",
                                         "GRS_IR_9wIR",
                                         "GRS_IS_17nw",
                                         "GRS_IS_17wIS",
                                         "GRS_LDL_54w",
                                         "GRS_LIPtotal_126nw",
                                         "GRS_T1D_29nw",
                                         "GRS_T1D_29w",
                                         "GRS_T2D_76nw",
                                         "GRS_T2D_76nw_IMP",
                                         "GRS_T2D_76nw_IMPtot",
                                         "GRS_TG_38w",
                                         "F_Incident_T2DM_022017",
                                         "F_Incident_T2DM_062016",
                                         "F_Incident_T2DM_072017",
                                         "F_Incident_T2DM_072018",
                                         "F_Incident_T2DM_092014",
                                         "F_Incident_T2DM_102018",
                                         "F_Incident_T2DM_2013",
                                         "F_FollowUpTime_Death_Feb2017",
                                         "F_FollowUpTime_T2DM_022017",
                                         "F_FollowUpTime_T2DM_062016",
                                         "F_FollowUpTime_T2DM_072017",
                                         "F_FollowUpTime_T2DM_072018",
                                         "F_FollowUpTime_T2DM_092014",
                                         "F_FollowUpTime_T2DM_102018",
                                         "F_FollowUpTime_T2DM_2013",
                                         "Death_Date_Dec2015",
                                         "Death_Date_Feb2017",
                                         "Death_Date_March2015"
))

metadata2$MOMID <- as.character(metadata2$MOMID)

metadata2$MOMID <- as.character(metadata2$MOMID)
d$METSIM_ID <- as.character(d$METSIM_ID)
d <- left_join(d, metadata2, by =c('METSIM_ID'='MOMID'))

d$X.1 <- NULL
d$X <- NULL
d$Sample_ID <- str_replace(d$SampleID, "[.]", "_")
write.csv(d, "FINAL_MICROBIOME_DATASET.csv")

######################## Single data
# Add new metabolomics data

targetedmetabolomics <- readxl::read_xlsx('METSIM plasma TMAO etc_053121.xlsx')
untargetedmetabolomics <- readxl::read_xlsx('METSIM_MetabolonData_Microbiome_converted_to_numeric.20200623.xlsx')


targetedmetabolomics$METSIM_ID <- as.character(targetedmetabolomics$METSIM_ID)
untargetedmetabolomics$metsim_id <- as.character(untargetedmetabolomics$metsim_id)
metadata$MOMID <- as.character(metadata$MOMID)
d <- inner_join(metadata, targetedmetabolomics, by = c("MOMID" = "METSIM_ID")) %>% left_join(untargetedmetabolomics, by = c("MOMID"="metsim_id"))
write.csv(d, file = 'FINAL_METABOLOMISC_DATASET.csv')
test <- left_join(metadata, untargetedmetabolomics, by=c('MOMID' = 'metsim_id')) %>% select(c("MOMID", "Res_date.x", "Res_date.y", 'F_Res_date','F2_Res_date', 'D0_0MonthBaselineResDate'))



############################# Check missing D0 and D18 #############################
datesdf <- readxl::read_xlsx("sample_dates.xlsx")
table(D18$MOMID_TP %in% NMR18full$MOMID_TP)
q <- which(D18$MOMID_TP %notin% NMR18full$MOMID_TP)
missings <- data.frame("MOMID" = D18$MOMID[q]) %>% left_join(datesdf, by = c("MOMID" = "METSIMID"))


table(D0$MOMID_TP %in% NMR0$MOMID_TP)
q <- which(D0$MOMID_TP %notin% NMR0$MOMID_TP)
missings2 <- data.frame("MOMID" = D0$MOMID[q]) %>% left_join(datesdf, by = c("MOMID" = "METSIMID"))

writexl::write_xlsx(list("T18" = missings,"T0" = missings2), "samples_missing_NMR_data.xlsx")

singlesample <- subset(datesdf, datesdf$`Total samples` ==1) %>% left_join(D01848microbiome, by = c('METSIMID'='METSIM_ID')) %>% left_join(NMR01848, by = 'MOMID_TP')
write.csv(singlesample, 'singlesampledataset.csv')


# Save to file 
write.csv(d, file="full_dataset.csv")


###########################################
# d <- read.csv('full_dataset.csv')
# s <- readxl::read_xlsx('sample_dates.xlsx', 1)
# d <- left_join(d, s, by = c('METSIM_ID' = "METSIMID"))
# d$`Total samples`[d$METSIM.ID == 2425] <- 2
# d$`Total samples`[d$METSIM.ID == 5414] <- 2
# save(d, file = 'full_dataset.rda')
# t <- subset(d, d$`Total samples`==2)
# t
# write.csv(t, file ='double_samples.csv')
# 
# table(is.na(t))
# 
# na_count <- data.frame(sapply(d, function(y) sum(length(which(is.na(y))))))
# na_count$variable <- rownames(na_count)
# na_count <- left_join(na_count, explainer, by = c('variable'='trait_name'))
# write.csv(na_count, 'NA_count_All.csv')
# 
# f <- readxl::read_xlsx('Metsim_metadata/Metadata_metsim_data_20181219.xlsx')
# explainer <- readxl::read_xlsx('Metsim_metadata/Metadata_metsim_data_20181219.xlsx',2)
# f2 <- inner_join(f, s, by = c("MOMID" = "METSIMID"))
# f2$`Total samples`[f2$MOMID == 2425] <- 2
# f2$`Total samples`[f2$MOMID == 5414] <- 2
# f3 <- subset(f2, f2$`Total samples` ==2)
# na_count3 <- data.frame(sapply(f3, function(y) sum(length(which(is.na(y))))))
# write.csv(na_count3, 'NA_count3_doubles.csv')
# write.csv(f3, 'double_samples_baseline_data.csv')
# f3 <- subset(f3, MOMID!= 2425 & MOMID !=5414)
# na_count3$var <- rownames(na_count3)
# na_count3 <- left_join(na_count3, explainer, by = c('var'='trait_name'))
# 
# 
# na_count <- data.frame(sapply(f2, function(y) sum(length(which(is.na(y))))))
# na_count$var <- rownames(na_count)
# na_count <- left_join(na_count, explainer, by = c('var'='trait_name'))
# write.csv(na_count, 'NA_count_all_baseline.csv')
