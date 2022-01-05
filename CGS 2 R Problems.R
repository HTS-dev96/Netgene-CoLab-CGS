## 1. Download the clinical data for TCGA-GBM

##  2. Make a new dataset from columns: 
## vital_status, gender, race, days_to_death, treatments_pharmaceutical_treatment_or_therapy, treatments_radiation_treatment_or_therapy

## 3. Use UMAP followed by PAM clustering, search for 2-3 clusters, do not include days to death in the UMAP
## optimize using gap-stat
## Rows 600-617 in the GBM dataset are empty
## Hint, don't forget to code the characters into numerical variables
## And if someone is alive then days to death is NA and need to be filled in
## I'd recommend 0 because it won't come up but you'll notice your mistake

## 4. Plot a KM survival curve based off the clusters

## 5. Plot a coxph with 5 variables of your choosing using the original dataset

##  6. Fit a cox model using glment on the dataset you made in #2