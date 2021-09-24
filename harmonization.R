rm(list=ls())
#Prep PEDI data 
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(neuroCombat)

#start with sheet containing demographics and tract-wise connectivity with one row per participant
dat <- read.csv("demographics_and_brain.csv")

original <- dat

#select tract data only
df <- t(dat[,c(40:3609)])
df <- as.data.frame(df)
colnames(df) <- c(dat$code)
df[,c(1:267)] <- sapply(df[,c(1:267)], as.character)
df[,c(1:267)] <- sapply(df[,c(1:267)], as.numeric)
df[is.na(df)] <- 0

#Map data to -inf to inf domain - prevent impossible values post-harmonisation
df[,c(1:267)] <- sapply(df[,c(1:267)], sqrt)

#Remove constant rows
df <- df[rowSums(is.na(df)) != ncol(df), ]
keep <- apply(df[1:267], 1, function(x) length(unique(x[!is.na(x)])) != 1)
df <- df[keep, ]

mat <- as.matrix(df)

#Set up harmonization model
BMI <- dat$BMI
Age <- dat$Age
Dx <- dat$Diagnosis

#batch should be numeric - this is number of different sites/scanners/studies
batch <- dat$Study


m <- model.matrix(~ Dx + Age + BMI)

combat.harmonized <- neuroCombat::neuroCombat(mat, batch=batch, mod=m)


## test whether there are still distributional differences in feature 1 between batches 1 and 2, post-combat
tracts <- rownames(mat)
har <- combat.harmonized$dat.combat
har <- t(har)
har <- as.data.frame(har)

#push harmonized and demographic data back together again
#select demographics from data file
dem <- dat[,c(1:38)]
dat <- merge(dem, har, by="ID")

#square your outputs from harmonization - reverse the transform
dat[,c(39:3521)] <- sapply(dat[,c(39:3521)], as.character)
dat[,c(39:3521)]<- sapply(dat[,c(39:3521)], as.numeric)
dat[,c(39:3521)] <- (dat[,c(39:3521)]^2)



#check whether harmonization did a good job of evening out variances 
D <- c()
for(i in 39:3523){
  n1 <- sample(1:5, 1)
  n2 <- sample(1:5, 1)
  res <- tryCatch({
  if(!(n1==n2)){
  fr <- data.frame(con=dat[[i]], batch=batch)
  v1 <- fr$con[batch==n1]
  v2 <- fr$con[batch==n2]
  p1 <- ks.test(v1, v2)$p.value
  fr <- data.frame(con=OG[[i]], batch=batch)
  v1 <- fr$con[batch==n1]
  v2 <- fr$con[batch==n2]
  p2 <- ks.test(v1, v2)$p.value
  if(p1 < 0.05 & p2 > 0.05){
    d <- ("harm fix")
  }  
  if(p1 > 0.05 & p2 < 0.05){
    d <- "harm prob"  
  }
  if(p1 > 0.05 & p2 > 0.05){
    d <- "no prob" 
  }
  if(p1 < 0.05 & p2 < 0.05){
    d <- "no fix"  
  }
  D <- c(D,d)  
  }}, error=function(e) e)
}

