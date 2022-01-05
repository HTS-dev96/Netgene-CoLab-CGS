if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install(c("TCGABiolinks"))

library(TCGAbiolinks)
library(SummarizedExperiment)
## Always convert from FKPM to TPM
## Use the function below
fpkmToTpm <- function(fpkm) {
  
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  
}

## Pulling RNA-seq data
## Need to runa query first
## Puts together all of your parameters and "finds" the data

query.rna.1 <- GDCquery(
  project = c("TCGA-ACC", 'TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC','TCGA-CHOL'),
  legacy = F,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = 'HTSeq - FPKM')


## Then you need to download it

GDCdownload(query.rna.1, method = 'api')

## Then you must save it
## save.filename is optional but it will save it outside of your R enviroment
## which is so beneficial that unless the file size prevents you, you should always do it

rna.1 <- GDCprepare(query.rna.1, save = T, 
                    save.filename = 'rna.1.rda')

## Then convert into a dataframe for easy manipulation

rna.1.df <- as.data.frame(assay(rna.1))

## For clinical data use
clinical <- GDCquery_clinic(project = "Project code", type = "clinical")


## Kaplan-meier curve

library(survival)
library(survminer)
library(survcomp)

## Generating fake dataset
total_time <- rnorm(200, mean = 200, sd = 50)
status <- floor(runif(200, min = 0, max = 2))
variable <- floor(runif(200, min = 0, max = 2))

data <- data.frame(status, variable, total_time)

## Preparing survival curve
## KM plots can only have one variable
## type refers to the type of censoring that happened with the data
## Remember that total_time refers to last known status of each person
## So either last known follow-up or death. Those times can be different 
## Remember KM plots can only handle catagorical variables not continuous 

km <- survfit(Surv(total_time, status, type = 'right') ~ variable, data = data)

## Now you need to plot it
## A lot of repition 
## Feel free to play around with risk.table and censor. 
## I normally plot this way because I think it is more professional looking
## Do NOT set censor to T when you have thousands of data points. R will crash!

ggsurvplot(km, data = data, pval = T,
                          risk.table = F, censor = F)

## Cox PH Model
## New data set
total_time <- rnorm(200, mean = 200, sd = 50)
status <- floor(runif(200, min = 0, max = 2))
variable <- floor(runif(200, min = 0, max = 2))
covariate1 <- floor(runif(200, min = 0, max = 2))
covariate2 <- rnorm(200, mean = 20, sd =2.5)

data <- data.frame(status, variable, total_time, covariate1, covariate2)

## Very similar set up, but now you include other variables

data.cox <- coxph(Surv(total_time, status, type = 'right') ~
                                  variable + covariate1 + covariate2,
                                data = data)
## Test the coxph assumption

cox.zph(data.cox)

## This gives you a report of the results
summary(data.cox)

## This is how you graph it, again can only graph categorical variables
## Continuous variables give you results in the summary
ggadjustedcurves(data.cox,
                 variable = "covariate1",
                 method = "average",
                 xlim = c(0, 400), 
                 title = "Example Cox Model", 
                 data = data)

## Glmnet package
## Survival curve

install.packages("glmnet")
library(glmnet)

## glment requires everything to be a matrix and all values to known
## Matrices must have the same type of variable inside of them so if a single
## cell has a character it will become a character/string matrix
## You must have a numerical matrix
## To check run str(your matrix name here)
## All factor variabel must be turned into individual variables. For example
## if you have a variable for ethnicity you need to change it into "caucasian", 
## "african", "asian" etc. 

## Using their fake dataset
data(CoxExample)
x <- CoxExample$x
y <- CoxExample$y
y[1:5, ]

## You use the function 'glment"
## The first input is variables
## The second input is the survial, look at the variable y to see how it is structured
## Third is type of regression

fit <- glmnet(x,
              y,
              family = 'cox')

## To visualize the fit
plot(fit)

## Now you want to cross-validate to find the optimal fit

fit.cv <- cv.glmnet(x,
                    y,
                    family = "cox",
                    type.measure = "C")
## Visualize it
plot(fit.cv)

## Now extract coefficients at optimal fit

fit.optimal <- coef(fit, s = fit.cv$lambda.min)

summary(fit.optimal)

## i is the order of the varible, j you can ignore, and x is the coefficient 