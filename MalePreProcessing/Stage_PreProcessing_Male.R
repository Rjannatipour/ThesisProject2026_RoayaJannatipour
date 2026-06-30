#Data Pre-Processing- Study SDY1276

library(dplyr)
library(gtsummary)
library(purrr)
library(janitor)
library(tidyr)
library(boot)
library(matrixStats)
library(caret)
library(readr)
library(tidyverse)

# DEFINE MYFILEPATH:

myfilepath <- "/Users/roayajannatipour/Desktop/ThesisProject2026_RoayaJannatipour"

# -----Create data objects-----
is2_metadata <- read_csv("rawdata/is2_metadata.csv") %>% 
  clean_names()
is2rawab <- read_csv("rawdata/is2rawab.csv") %>% 
  clean_names()

is2norm = get(load("/Users/roayajannatipour/Desktop/ThesisProject2026_RoayaJannatipour/rawdata/IS2_all_norm_expression(1).RData"))


#----Data Cleaning -----
#Select participants from the 'SDY1276' study located in 'is2_metadata'
#Creating data intersection
study_data = is2_metadata %>% 
  select(participant_id, study_accession) %>% 
  filter(study_accession == 'SDY1276') %>% 
  select(participant_id)%>% 
  unique()

study_ab = is2rawab %>% 
  select(participant_id) 

israwab_study = intersect(study_data, study_ab) # Select participant ids that are in both in the metadata and have recorded antibody response under 'value_preferred'

nrow(israwab_study) #218 participant_ids in this intersection


is2rawab_study = israwab_study %>% 
  left_join(is2rawab) %>% 
  drop_na() #217 after left joining participant_ids with dataframe containing antibody response measures


#Filter for participants with both Visit 0 and Visit 28 recorded

pre_is2rawab_study = is2rawab_study %>% 
  filter(visit == 0)

pid_pre = pre_is2rawab_study %>% 
  pull(participant_id)

post_is2rawab_study = is2rawab_study %>% 
  filter(visit == 28)

pid_post = post_is2rawab_study %>% 
  pull(participant_id)

pids_both <- intersect(pid_pre, pid_post)

pre_is2rawab_study = pre_is2rawab_study %>% 
  filter(participant_id %in% pids_both) %>% 
  arrange(participant_id, virus)

post_is2rawab_study = post_is2rawab_study %>% 
  filter(participant_id %in% pids_both) %>% 
  arrange(participant_id, virus)

is2rawab_study_wide = is2rawab_study %>% 
  filter(visit %in% c(0,28),
         participant_id %in% pids_both) %>% 
  pivot_wider(id_cols = c(participant_id, virus),
              names_from = visit,
              values_from = value_preferred,
              names_prefix = "visit_"
  ) %>% 
  filter(!is.na(visit_0 | visit_28))

is2rawab_study_wide  %>% 
  summarise(n_unique = n_distinct(`participant_id`))

rm(list = setdiff(ls(), c("is2rawab_study_wide", "pids_both", "is2_metadata", "is2rawab", "is2norm")))

#214 unique participants, who have both visit 0 and visit 28, in this study. 
#Each participant has visit 0 and visit 28 from 3 different strands of virus
#Next step is to determine the max FC per participant id

#-------Creating Base Model------- 
#Join metadata to is2rawab_study_wide
basemodelmeta = is2_metadata %>% select(participant_id, age_imputed) %>% distinct()
basemodelab = is2rawab%>%  select(participant_id, gender) %>% distinct()
BaseModel = is2rawab_study_wide %>%  left_join(basemodelmeta, by= "participant_id")  %>%  left_join(basemodelab, by = "participant_id") %>% na.omit()

BaseModel  %>% distinct(participant_id, .keep_all=TRUE) %>% tabyl(gender) %>%  adorn_pct_formatting() %>% flextable::flextable()
#106 females and 108 males

BaseModel <- BaseModel %>% filter(gender == 'Male')
#Only male participants in our basemodel (108) #Calculating max FC for each participant

BaseModel = BaseModel %>% group_by(participant_id) %>% mutate(maxFC = max(visit_28-visit_0)) %>% slice_max(visit_28 - visit_0, n = 1, with_ties = FALSE) %>% ungroup()

BaseModel %>% distinct(participant_id) %>% nrow() #108 participants (males) in Base Model

rm(basemodelmeta, basemodelab)

#------Data Cleaning for Univariate Model------

valid_ids <- is2norm %>%
  filter(study_time_collected %in% c(0,1)) %>%
  group_by(participant_id) %>%
  summarise(
    n0 = sum(study_time_collected == 0),
    n1 = sum(study_time_collected == 1),
    .groups = "drop"
  ) %>%
  filter(n0 == 1, n1 == 1) %>%
  pull(participant_id)


is2norm01 = is2norm %>% 
  filter(participant_id %in% valid_ids,
    study_time_collected %in% c(0,1)) %>% 
  group_by(participant_id)%>%
  summarise(
    across(
      a1cf:zzz3,
      ~ .[study_time_collected == 1] - .[study_time_collected == 0]
    ),
    .groups = "drop"
  )

is2norm01 %>% distinct(participant_id) %>%nrow()
#477 distinct participant IDs with genetic expression data at Day 0 and Day 1

#Merged Data Frame

merge<- BaseModel %>% select(participant_id,age_imputed,maxFC)%>% left_join(is2norm01, by = "participant_id" ) %>% drop_na() 

merge  %>% distinct(participant_id, .keep_all=TRUE) %>% nrow()

#105 males in merged data

#----Intersect BaseModel and IntModel Participants to discard 3 participants w/o genetic expression data---

FinalParticipants <- merge %>% select(participant_id)

BaseModel <- BaseModel %>% right_join(FinalParticipants, by = 'participant_id') %>% arrange(participant_id) #105 males in Base Model

saveRDS(BaseModel, "MalePreProcessing/BaseModel.rds") 
saveRDS(is2norm01, "rawdata/is2norm01.rds")
saveRDS(merge,"rawdata/mergemale.rds")

#---Creating RISE data objects----

#Filter genetic expresion data to only collect antibody values at time 0

is2normrise = is2norm %>% 
  filter(participant_id %in% valid_ids) 

is2norm_studytime0_male = is2normrise %>% 
  filter(study_time_collected == 0)

is2norm_studytime0_male %>% distinct(participant_id) %>%nrow()
##477 distinct participant IDs with study time collected 0

is2norm_studytime1_male = is2normrise %>% 
  filter(study_time_collected == 1)

is2norm_studytime1_male %>% distinct(participant_id) %>%nrow()
##477 distinct participant IDs with study time collected 1


#Create dataframes to merge our male participants with their corresponding genetic data at Day 0

merge_studytime0_male <- BaseModel %>% 
  select(participant_id)%>% 
  left_join(is2norm_studytime0_male, by = "participant_id" ) %>%
  drop_na() %>% 
  arrange(participant_id) %>% 
  as.data.frame()

row.names(merge_studytime0_male) <- merge_studytime0_male$participant_id

merge_studytime0_male <-  merge_studytime0_male %>% select(-study_time_collected)

merge_studytime0_male  %>% distinct(participant_id, .keep_all=TRUE) %>% nrow()

#105 males in merged data, with gene expression data for visit 0

#Create dataframes to merge our male participants with their corresponding genetic data at Day 1

merge_studytime1_male <- BaseModel %>% 
  select(participant_id)%>% 
  left_join(is2norm_studytime1_male, by = "participant_id" ) %>% 
  drop_na() %>% 
  arrange(participant_id) %>%
  as.data.frame()

row.names(merge_studytime1_male) <- merge_studytime1_male$participant_id

merge_studytime1_male <-  merge_studytime1_male %>% select(-study_time_collected)

merge_studytime1_male  %>% distinct(participant_id, .keep_all=TRUE) %>% nrow()

identical(merge_studytime1_male$participant_id, merge_studytime0_male$participant_id)

#Identical lists of male participants

#105 males in merged data, with gene expression data for visit 1

# Remove participant_id column from both merged dataframes

merge_studytime0_male <-  merge_studytime0_male %>% select(-participant_id)

merge_studytime1_male <-  merge_studytime1_male %>% select(-participant_id)


#Assign RISE specific variables to our data objects

yonemale <- BaseModel$visit_28

yzeromale <- BaseModel$visit_0

szeromale <- merge_studytime0_male

sonemale <- merge_studytime1_male

#Save RISE objects

saveRDS(yonemale, "MalePreProcessing/yonemale.rds")
saveRDS(yzeromale, "MalePreProcessing/yzeromale.rds")
saveRDS(szeromale, "MalePreProcessing/szeromale.rds")
saveRDS(sonemale, "MalePreProcessing/sonemale.rds")


saveRDS(merge_studytime0_male, "MalePreProcessing/merge_studytime0_male.rds")
saveRDS(merge_studytime1_male, "MalePreProcessing/merge_studytime1_male.rds")


#Mixed Model

BaseModel_Male_Long <- BaseModel %>%
  pivot_longer(
    cols = c(visit_0, visit_28),
    names_to = "visit",
    values_to = "value"
  ) %>% 
  mutate(Visit = gsub("visit_", "", visit)) %>% 
  mutate(Vaccinated = if_else(Visit >0, 1, 0), Antibody = value) %>% select(-value) %>%
  mutate(Age = age_imputed) %>% select(-age_imputed)

BaseModel_Male_Long  <- BaseModel_Male_Long  %>%
  mutate(Vaccinated = as.factor(Vaccinated))

mixedmale <- lmer(Antibody ~ Age + Vaccinated+ (1|participant_id), BaseModel_Male_Long, REML = 0)

mixedmaleint <- lmer(Antibody ~ Age * Vaccinated + (1 | participant_id),data = BaseModel_Male_Long, REML = 0)

#Mixed model residuals

BaseModel_Male_Long$condres <- resid(mixedmale)


# RISE Prep

yonemalemix <- BaseModel_Male_Long %>% filter(Visit == 28) %>% pull(condres) 

yzeromalemix <- BaseModel_Male_Long %>% filter(Visit == 0) %>% pull(condres)

saveRDS(BaseModel_Male_Long, "MalePreProcessing/BaseModel_Male_Long.rds")
saveRDS(yonemalemix, "MalePreProcessing/yonemalemix.rds")
saveRDS(yzeromalemix, "MalePreProcessing/yzeromalemix.rds")
