#FLOWERING TIME####

setwd('/EpiQTL/Flowering_time')

library(qtl)
dat <- readRDS('Population_analysis_Flowering_Time.RDS')
df1 <- NULL
for ( i in 1:length(dat)){
  hvdrr <- dat[[i]]
  hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Flowering')))
  hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  
  #replace 3 with 2
  hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  class(hvdrr)[1] <- 'riself'
  summary(hvdrr)
  hvdrrprob2 <- calc.genoprob(hvdrr, step = 2.5, err = 0.01, map.function = "haldane")
  hvdrrout2 <- scantwo(hvdrrprob2, verbose=FALSE)
  summary(hvdrrout2)
  df.s2sum <- data.frame(summary(hvdrrout2))
  df.s2sum$populaiton <- names(dat[i])
  df1 <- rbind(df1, df.s2sum)
}
write.csv(df1, 'epiQTL.FT2.5.csv', row.names = F)

#filter for epiQTL with LOD > 3

setwd('/EpiQTL/Flowering_time')
#setwd('/home/shrestha/mounts/project/shrestha/marvin/seedsizeQTL/experiments/qtlmapping/')
library(qtl)
dat <- readRDS('Population_analysis_Flowering_Time.RDS')
df1 <- read.csv('epiqtl_FT_lod3.csv')
df1.ls <- levels(df1$populaiton)

#i <- 4

for (i in 1:3){
  hvdrr <- dat[[df1.ls[i]]]
  print(summary(hvdrr))
  hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Flowering')))
  
  hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  
  #replace 3 with 2
  hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  
  class(hvdrr)[1] <- 'riself'
  
  hvdrrprob2 <- calc.genoprob(hvdrr, step = 2.5, err = 0.01, map.function = "haldane")
  set.seed(85842100)
  out.perm2A <- scantwo(hvdrrprob2, n.perm = 500, clean.output=TRUE)
  save(out.perm2A, file = paste('/Epistasis_with_permutations/','','operm.',df1.ls[i],'.RData', sep = ''))
}


#PLANT HEIGHT####

setwd('/EpiQTL/Plant_height')
library(qtl)
dat <- readRDS('Population_analysis_Plant_Height.RDS')
df1 <- NULL
for ( i in 1:length(dat)){
  hvdrr <- dat[[i]]
  hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Plant_Height')))
  hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  
  #replace 3 with 2
  hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  class(hvdrr)[1] <- 'riself'
  summary(hvdrr)
  hvdrrprob2 <- calc.genoprob(hvdrr, step = 2.5, err = 0.01, map.function = "haldane")
  hvdrrout2 <- scantwo(hvdrrprob2, verbose=FALSE)
  summary(hvdrrout2)
  df.s2sum <- data.frame(summary(hvdrrout2))
  df.s2sum$populaiton <- names(dat[i])
  df1 <- rbind(df1, df.s2sum)
}
write.csv(df1, 'epiQTL.PH2.5.csv', row.names = F)

#filter for epiQTL with LOD > 3


setwd('/EpiQTL/Plant_height/')
#setwd('/home/shrestha/mounts/project/shrestha/marvin/seedsizeQTL/experiments/qtlmapping/')
library(qtl)
dat <- readRDS('Population_analysis_Plant_Height.RDS')
df1 <- read.csv('epiqtl_PH_lod3.csv')
df1.ls <- levels(df1$populaiton)

for (i in 1:4){
  hvdrr <- dat[[df1.ls[i]]]
  print(summary(hvdrr))
  hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Plant_Height')))
  
  hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
    ifelse(x == 2, NA, x)
  })
  
  #replace 3 with 2
  hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
    ifelse(x == 3, 2, x)
  })
  
  class(hvdrr)[1] <- 'riself'
  
  hvdrrprob2 <- calc.genoprob(hvdrr, step = 2.5, err = 0.01, map.function = "haldane")
  set.seed(85842100)
  out.perm2A <- scantwo(hvdrrprob2, n.perm = 500, clean.output=TRUE)
  save(out.perm2A, file = paste('/Epistasis_with_permutations/','','operm.',df1.ls[i],'.RData', sep = ''))
}
