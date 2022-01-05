## 1. Download the clinical data for TCGA-GBM
library(TCGAbiolinks)

clinical <- GDCquery_clinic(project = "TCGA-GBM", type = "clinical")

##  2. Make a new dataset from columns: 
## vital_status, gender, race, days_to_death, treatments_pharmaceutical_treatment_or_therapy, treatments_radiation_treatment_or_therapy
library(tidyverse)

clinical.1 <- clinical %>%
  select(vital_status,
         gender,
         race,
         days_to_death,
         treatments_pharmaceutical_treatment_or_therapy,
         treatments_radiation_treatment_or_therapy)

## 3. Use UMAP followed by PAM clustering, search for 2-3 clusters, do not include days to death in the UMAP
## optimize using gap-stat
## Rows 600-617 in the GBM dataset are empty
## Hint, don't forget to code the characters into numerical variables
## And if someone is alive then days to death is NA and need to be filled in
## I'd recommend 0 because it won't come up but you'll notice your mistake

library(umap)
library(cluster)
library(factoextra)
## Dead is 1, alive is 0
## Female is 1, male is 0
## Asian is 3, black or african american is 2, not reported is 1, white is 0
## Pharmaceutical treatment Yes is 3, no is 2, not reported 1, NA is 0
## Rad treatment is the same as above

clinical.1$vital_status <- ifelse(clinical.1$vital_status == "Dead", 1, 0)
clinical.1$gender <- ifelse(clinical.1$gender == 'male', 0, 1)
clinical.1$race <- ifelse(clinical.1$race == "asian", 3,
                          ifelse(clinical.1$race == "black or african american", 2,
                                 ifelse(clinical.1$race == 'not reported', 1, 0)))

clinical.1$treatments_pharmaceutical_treatment_or_therapy <- ifelse(clinical.1$treatments_pharmaceutical_treatment_or_therapy == "yes", 3,
                                                                    ifelse(clinical.1$treatments_pharmaceutical_treatment_or_therapy == "no", 2,
                                                                           ifelse(clinical.1$treatments_pharmaceutical_treatment_or_therapy == 'not reported', 1, 0)))

clinical.1$treatments_radiation_treatment_or_therapy <- ifelse(clinical.1$treatments_radiation_treatment_or_therapy == "yes", 3,
                                                               ifelse(clinical.1$treatments_radiation_treatment_or_therapy == "no", 2,
                                                                      ifelse(clinical.1$treatments_radiation_treatment_or_therapy == 'not reported', 1, 0)))

clinical.1$days_to_death <- ifelse(is.na(clinical.1$days_to_death), 0, clinical.1$days_to_death)

clinical.1 <- clinical.1[1:599,]

umap <- umap(clinical.1[,c(1:3, 5, 6)], preserve.seed = T)

umap.df <- as.data.frame(umap[['layout']])

fviz_nbclust(umap.df,
             FUNcluster = cluster::pam,
             method = 'gap_stat',
             k.max = 3,
             print.summary = T,
             nboot = 50)

clinical.clusters <- pam(clinical.1[,c(1:3, 5, 6)], 3)

clinical.clus <- clinical.clusters[['clustering']]

clinical.1$clusters <- clinical.clus

## 4. Plot a KM survival curve based off the clusters
library(survival)
library(survminer)
library(survcomp)

km <- survfit(Surv(days_to_death, vital_status, type = 'right') ~ clusters, data = clinical.1)

ggsurvplot(km, data = clinical.1, pval = T,
           risk.table = F, censor = F)

## 5. Print a summary of a coxph with 5 variables of your choosing using the original dataset
## I'm going to use the modified dataset but it all works the same

cox <- coxph(Surv(days_to_death, vital_status, type = 'right') ~ gender + race + 
                 treatments_pharmaceutical_treatment_or_therapy, data = clinical.1)

summary(cox)

##  6. Fit a cox model using glment on the dataset you made in #2

## It's late I hope I did this later lol
