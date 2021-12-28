## Introduction to R

## A hashtag tells R to ignore what comes after it. You can use one or two
## to leave notes for yourself. This is your lab notebook 

## To run lines of code press command + enter at the same time.

## R is object based so you and one of the most useful objects are dataframes
## Setting the seed ensures that we both get the same results if we do something
## that invovles random number generation
set.seed(1014)

## To save an object use <-
## We will use the example dataset 'mtcars'
data("mtcars")
## To look at your df you can 
View(mtcars)
## But if it is large you might want to use 
head(mtcars)
tail(mtcars)

## For data wrangling the best package to use is 'tidyverse'
## The data.table package can run code in parallel (which just means it runs faster)
## and shares many of the same functions as tidyverse

## To install a package from the base R respository
install.packages('tidyverse')
install.packages('data.table')
## Now you need to tell R you want to use that package now
library(data.table)
library(tidyverse)

## Now that it is installed you don't have to install it again
## the tidyverse has many different packages in it
## many cheatsheets are available online which I included in the email
## one of the most important operators to learn is '%>%'
## It is called piping and allows you to easily nest functions which allows for
## easy data wrangling

## For example what if you want to pull mpg, and wt data on cars with 4 gears?
## == is used for equals, != is not equal to 
## select chooses columns
df1 <- mtcars %>%
  select(mpg, wt, gear) %>%
  filter(gear == 4)

## With many variables you can make life a little easier with 
starts_with()
ends_with()

## To have R recognize something as characters it must be between ' or "
df2 <- mtcars %>%
  select(starts_with('TCGA'))

## You can also create new columns with the function mutate

df3 <- mtcars %>%
  mutate(V1 = mpg/wt)

## Another very important feature is being able to filter by checking for values
## within a list such as gene names
## c() is how to tell R that multiple things will be in there (numbers or characters)
genes <- c('MYC', 'LDHA', 'LDHB')

## Don't run, just for example


genes_of_interest <- rna_seq %>%
  filter(hgnc_id %in% genes)

## Then you get a list of your genes of interest and their counts from an RNA-seq

## An amazing package in R is TCGAbiolinks
## It comes from a biology focused package repository called Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
library('TCGAbiolinks')

## It can allow you to pull data from any TCGA project and more which allows you
## to compare across tumors and have a large sample

## Here is how to read up on it
browseVignettes("TCGAbiolinks")

## Mainly you'll need 3 functions
## To learn more about these functions and any function type a '?' in front
## before running it

?GDCquery

GDCquery()

GDCdownload()

GDCprepare()

## There are many catagories for these functions
## Then to pull the assay data you'll need to use the package summarized experiments
BiocManager::install('SummarizedExperiment')

## If you use this function then assay() will pull out the RNA seq data
GDCprepare_output <- GDCprepare(query_RNA_seq)
df <- assay(GDCprepare_output)

## How to do clustering
## There are many types of clustering. They are well explained with this website
## https://towardsdatascience.com/the-5-clustering-algorithms-data-scientists-need-to-know-a36d136ef68#:~:text=Given%20a%20set%20of%20data,dissimilar%20properties%20and%2For%20features.

## General take aways
## K-meoids is an improvement on K-means because it allows for a variety of shapes
## CLARA clustering is optimal for very large data sets. Samples >> 2,000
## Hierarchical clustering can be really useful if there is an underlying hierarchy
## you are trying to recover. 
## Lastly, hierarchical clustering will find the optimal number of clusters other methods
## Will require testing to find


## One of the best packages I've found is cluster
install.packages('clutser')
library(cluster)
library(factoextra)
## we will start with partitional clustering (K-means etc)
## first find the optimal number of clusters

fviz_nbclust(ras.matrix.transpose,
             FUNcluster = cluster::## clustering method you want,
             method = ## 3 choices here, hard to tell which is best I like gap_stat,
             k.max = ## maximum number of clusters to test,
             print.summary = T,
             nboot = ## Only need for 'gap_stat' method
               )

## as an example
data("iris")
## need to remove any characters
## PAM is also known as K-meoids
fviz_nbclust(iris[,1:4],
             FUNcluster = cluster::pam,
             method = 'gap_stat',
             k.max = 10,
             print.summary = T,
             nboot = 50)
## Optimal number of clusters is 4

## Now do the clustering
## CLARA requires you to pick the size of samples because it makes mini clusters before expanding
iris.clusters <- pam(iris[,1:4], 4)
## Lets look!
fviz_cluster(iris.clusters, repel = T, geom = 'point', pointsize = .75)
## It automatically does a PCA and graphs the first two principal components
## While it is nice to have principle components it is not the best way to visuaulize
## UMAP (the better version) is covered later

## Pulling out the clusters

iris.clusters.df <- iris.clusters[['clustering']]
iris$clusters <- iris.clusters.df

## Visualization
## Many visualization methods, t-sne, PCA, and UMAP
## UMAP is becoming very popular and can actually improve clustering accuracy
## Others do not increase clustering accuary
## UMAP works well with linear and non-linear relationships and global/local relationships

install.package('umap')
library(umap)

## Only numeric variables can be used

iris.umap <- umap(iris[,1:4], preserve.seed = T)
## Pull the UMAP results and turn it into a dataframe
## It needs to be a dataframe to make graphing it easier
iris.umap.df <- as.data.frame(iris.umap[['layout']])
## Give the columns names
colnames(iris.umap.df) <- c("UMAP_1", "UMAP_2")

ggplot(iris.umap.df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point() +
  ggtitle("UMAP of Iris Data")

## Now lets color the points based on our clusters
## we need to combine the datasets to do this 
## because we haven't altered the datasets yet they are still in the same order
iris.umap.clusters.df <- cbind(iris, iris.umap.df)

ggplot(iris.umap.clusters.df, aes(x = UMAP_1, y = UMAP_2, color = clusters)) +
  geom_point()
## It is hard to see the blue but it is the default 
## to change it 
install.packages('RColorBrewer')
library(RColorBrewer)
ggplot(iris.umap.clusters.df, aes(x = UMAP_1, y = UMAP_2, color = clusters)) +
  geom_point() +
  scale_color_brewer(palette="spectral")

## Survival Analysis 

## Kaplan Meier Curve
## Always a good idea to do first because there are no assumptions unlike the Cox
## PH model. 
survival <- c('survival', 'survivalmodels', 'survcomp', "survMisc")
install.package(survival)
library(survival)

## KM





