#FLOWERING TIME##################
setwd("/QTLanalysis/Flowering_time/")

#HvDRRO2#####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR02 <- Populations$HvDRR02

HvDRR02$geno$'1'$data<- apply(HvDRR02$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'2'$data<- apply(HvDRR02$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'3'$data<- apply(HvDRR02$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'4'$data<- apply(HvDRR02$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'5'$data<- apply(HvDRR02$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'6'$data<- apply(HvDRR02$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'7'$data<- apply(HvDRR02$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
#to eliminate heterozygous
HvDRR02$geno$'1'$data<- apply(HvDRR02$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'2'$data<- apply(HvDRR02$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'3'$data<- apply(HvDRR02$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'4'$data<- apply(HvDRR02$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'5'$data<- apply(HvDRR02$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'6'$data<- apply(HvDRR02$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'7'$data<- apply(HvDRR02$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR02)[1] <- "riself"



summary(HvDRR02)
plot(HvDRR02)
HvDRR02 <- subset(HvDRR02,ind=c(1:nrow(HvDRR02$pheno))[!is.na(HvDRR02$pheno)])

HvDRR02per <- calc.genoprob(HvDRR02, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR02out <- scanone(HvDRR02per, method = "hk")
plot(HvDRR02out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR02permout <- scanone(HvDRR02per, method = "hk", n.perm=4000)
summary(HvDRR02out, perms = HvDRR02permout, alpha = 0.05, pvalues = TRUE)

HvDRR02qtl1 <- makeqtl(HvDRR02per, chr = c(2,5), pos = c(61.1, 202.9), what = "prob")
#to create the 1st QTL object
HvDRR02outc5 <- addqtl(HvDRR02per, qtl = HvDRR02qtl1, method = "hk")
summary(HvDRR02outc5, perms = HvDRR02permout, alpha = 0.05, pvalues = TRUE)
#there were no LOD peaks above the threshold, so two QTLs

HvDRR02rqtl <- refineqtl(HvDRR02per, qtl = HvDRR02qtl1, method = "hk", verbose = FALSE)
summary(HvDRR02rqtl)
summary(fitqtl(HvDRR02per, qtl = HvDRR02rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR02per, qtl = HvDRR02rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2))


find.marker(HvDRR02per, 2, 56.3)
effectplot(HvDRR02per, mname1 = "JHI-Hv50k-2016-73780")
summary(HvDRR02per, mname1 = "JHI-Hv50k-2016-73780")
plotPXG(HvDRR02per, "JHI-Hv50k-2016-73780")

find.marker(HvDRR02per, 5, 205.3)
effectplot(HvDRR02per, mname1 = "JHI-Hv50k-2016-336770")
summary(HvDRR02per, mname1 = "JHI-Hv50k-2016-336770")
plotPXG(HvDRR02per, "JHI-Hv50k-2016-336770")


HvDRR02CI2 <- lodint(HvDRR02rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR02CI5 <- lodint(HvDRR02rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)

pt(-6.633, 133, lower.tail = FALSE)
#p<0.01
pt(14.925, 133, lower.tail = FALSE)
#p<0.01

#HvDRR03#####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR03 <- Populations$HvDRR03

HvDRR03$geno$'1'$data<- apply(HvDRR03$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'2'$data<- apply(HvDRR03$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'3'$data<- apply(HvDRR03$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'4'$data<- apply(HvDRR03$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'5'$data<- apply(HvDRR03$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'6'$data<- apply(HvDRR03$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'7'$data<- apply(HvDRR03$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'1'$data<- apply(HvDRR03$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'2'$data<- apply(HvDRR03$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'3'$data<- apply(HvDRR03$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'4'$data<- apply(HvDRR03$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'5'$data<- apply(HvDRR03$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'6'$data<- apply(HvDRR03$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'7'$data<- apply(HvDRR03$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR03)[1] <- "riself"



summary(HvDRR03)
plot(HvDRR03)
HvDRR03 <- subset(HvDRR03,ind=c(1:nrow(HvDRR03$pheno))[!is.na(HvDRR03$pheno)])


HvDRR03per <- calc.genoprob(HvDRR03, step = 1,  error.prob = 0.01, map.function = "haldane")
HvDRR03out <- scanone(HvDRR03per, method = "hk")
plot(HvDRR03out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR03permout <- scanone(HvDRR03per, method = "hk", n.perm = 4000)
summary(HvDRR03out, perms = HvDRR03permout, alpha = 0.05, pvalues = TRUE)

HvDRR03qtl1 <- makeqtl(HvDRR03per, chr = 5, pos = 242, what = "prob")
HvDRR03outc5 <- addqtl(HvDRR03per, qtl = HvDRR03qtl1, method = "hk")
summary(HvDRR03outc5, perms = HvDRR03permout, alpha = 0.05, pvalues = TRUE)
HvDRR03qtl2 <- makeqtl(HvDRR03per, chr =  c(5,4), pos = c(242,73.7), what = "prob")
HvDRR03outc54 <- addqtl(HvDRR03per, qtl = HvDRR03qtl2, method = "hk")
summary(HvDRR03outc54, perms = HvDRR03permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR03out, HvDRR03outc54, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)

find.marker(HvDRR03per, 5, 242)
effectplot(HvDRR03per, mname1 = "JHI-Hv50k-2016-335902")
summary(HvDRR03per, mname1 = "JHI-Hv50k-2016-335902")
plotPXG(HvDRR03per, "JHI-Hv50k-2016-335902")

find.marker(HvDRR03per, 4, 73.7)
effectplot(HvDRR03per, mname1 = "JHI-Hv50k-2016-234502")
summary(HvDRR03per, mname1 = "JHI-Hv50k-2016-234502")
plotPXG(HvDRR03per, "JHI-Hv50k-2016-234502")

HvDRR03rqtl <- refineqtl(HvDRR03per, qtl = HvDRR03qtl2, method = "hk", verbose = FALSE)
summary(HvDRR03rqtl)
summary(fitqtl(HvDRR03per, qtl = HvDRR03rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR03per, qtl = HvDRR03rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1+Q2))

HvDRR03CI5 <- lodint(HvDRR03rqtl, chr = 5, qtl.index = 1)
HvDRR03CI4 <- lodint(HvDRR03rqtl, chr = 4, qtl.index = 2)

pt(10.039, 66, lower.tail = FALSE)
#p<0.01
pt(-4.244, 66, lower.tail = FALSE)
#p<0.01


#HvDRR07####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR07 <- Populations$HvDRR07

HvDRR07$geno$'1'$data<- apply(HvDRR07$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'2'$data<- apply(HvDRR07$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'3'$data<- apply(HvDRR07$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'4'$data<- apply(HvDRR07$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'5'$data<- apply(HvDRR07$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'6'$data<- apply(HvDRR07$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'7'$data<- apply(HvDRR07$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'1'$data<- apply(HvDRR07$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'2'$data<- apply(HvDRR07$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'3'$data<- apply(HvDRR07$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'4'$data<- apply(HvDRR07$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'5'$data<- apply(HvDRR07$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'6'$data<- apply(HvDRR07$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'7'$data<- apply(HvDRR07$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR07)[1] <- "riself"

summary(HvDRR07)
plot(HvDRR07)
HvDRR07 <- subset(HvDRR07,ind=c(1:nrow(HvDRR07$pheno))[!is.na(HvDRR07$pheno)])

HvDRR07per <- calc.genoprob(HvDRR07, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR07out <- scanone(HvDRR07per, method = "hk")
plot(HvDRR07out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR07permout <- scanone(HvDRR07per, method = "hk", n.perm = 4000)
summary(HvDRR07out, perms = HvDRR07permout, alpha = 0.05, pvalues = TRUE)
HvDRR07qtl1 <- makeqtl(HvDRR07per, chr = c(2, 3), pos = c(44.3, 154), what = "prob")
HvDRR07outc2 <- addqtl(HvDRR07per, qtl = HvDRR07qtl1, method = "hk")
summary(HvDRR07outc2, perms = HvDRR07permout, alpha = 0.05, pvalues = TRUE)

HvDRR07rqtl <- refineqtl(HvDRR07per, qtl = HvDRR07qtl1, method = "hk", verbose = FALSE)
summary(HvDRR07rqtl)
summary(fitqtl(HvDRR07per, qtl = HvDRR07rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR07per, 2, 43.8)
effectplot(HvDRR07per, mname1 = "BK_14")
summary(HvDRR07per, mname1 = "BK_14")
plotPXG(HvDRR07per, "BK_14")

find.marker(HvDRR07per, 3, 155)
effectplot(HvDRR07per, mname1 = "SCRI_RS_150944")
summary(HvDRR07per, mname1 = "SCRI_RS_150944")
plotPXG(HvDRR07per, "SCRI_RS_150944")

summary(fitqtl(HvDRR07per, qtl = HvDRR07rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 +Q2))

HvDRR07CI2 <- lodint(HvDRR07rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR07CI3 <- lodint(HvDRR07rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)

pt(-10.586, 95, lower.tail = FALSE)
#p<0.01
pt(-5.931, 95, lower.tail = FALSE)
#p<0.05


#HvDRR08#####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR08 <- Populations$HvDRR08

HvDRR08$geno$'1'$data<- apply(HvDRR08$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'2'$data<- apply(HvDRR08$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'3'$data<- apply(HvDRR08$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'4'$data<- apply(HvDRR08$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'5'$data<- apply(HvDRR08$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'6'$data<- apply(HvDRR08$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'7'$data<- apply(HvDRR08$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'1'$data<- apply(HvDRR08$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'2'$data<- apply(HvDRR08$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'3'$data<- apply(HvDRR08$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'4'$data<- apply(HvDRR08$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'5'$data<- apply(HvDRR08$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'6'$data<- apply(HvDRR08$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'7'$data<- apply(HvDRR08$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR08)[1] <- "riself"

summary(HvDRR08)
plot(HvDRR08)
HvDRR08 <- subset(HvDRR08,ind=c(1:nrow(HvDRR08$pheno))[!is.na(HvDRR08$pheno)])


HvDRR08per <- calc.genoprob(HvDRR08, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR08out <- scanone(HvDRR08per, method = "hk")
plot(HvDRR08out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR08permout <- scanone(HvDRR08per, method = "hk", n.perm = 4000)
summary(HvDRR08out, perms = HvDRR08permout, alpha = 0.05, pvalues = TRUE)
HvDRR08qtl1 <- makeqtl(HvDRR08per, chr = 2, pos = 58.4, what = "prob")
HvDRR08outc2 <- addqtl(HvDRR08per, qtl = HvDRR08qtl1, method = "hk")
summary(HvDRR08outc2, perms = HvDRR08permout, alpha = 0.05, pvalues = TRUE)
HvDRR08qtl2 <- makeqtl(HvDRR08per, chr =  c(2, 7), pos = c(58.4, 343), what = "prob")
HvDRR08outc27 <- addqtl(HvDRR08per, qtl = HvDRR08qtl2, method = "hk")
summary(HvDRR08outc27, perms = HvDRR08permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR08out, HvDRR08outc27, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)

HvDRR08rqtl <- refineqtl(HvDRR08per, qtl = HvDRR08qtl2, method = "hk", verbose = FALSE)
summary(HvDRR08rqtl)
summary(fitqtl(HvDRR08per, qtl = HvDRR08rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR08per, qtl = HvDRR08rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1+Q2))

HvDRR08CI2 <- lodint(HvDRR08rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR08CI7 <- lodint(HvDRR08rqtl, chr = 7, qtl.index = 2, expandtomarkers = TRUE)

find.marker(HvDRR08per, 2, 58.4)
effectplot(HvDRR08per, mname1 = "JHI-Hv50k-2016-73500")
summary(HvDRR08per, mname1 = "JHI-Hv50k-2016-73500")
plotPXG(HvDRR08per, "JHI-Hv50k-2016-73500")

find.marker(HvDRR08per, 7, 342.5)
effectplot(HvDRR08per, mname1 = "SCRI_RS_175551")
summary(HvDRR08per, mname1 = "SCRI_RS_175551")
plotPXG(HvDRR08per, "SCRI_RS_175551")

pt(10.555, 92, lower.tail = FALSE)
#p<0.01
pt(-4.11, 92, lower.tail = FALSE)
#p<0.01


#HvDRR09####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR09 <- Populations$HvDRR09

HvDRR09$geno$'1'$data<- apply(HvDRR09$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'2'$data<- apply(HvDRR09$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'3'$data<- apply(HvDRR09$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'4'$data<- apply(HvDRR09$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'5'$data<- apply(HvDRR09$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'6'$data<- apply(HvDRR09$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'7'$data<- apply(HvDRR09$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
#to eliminate heterozygous
HvDRR09$geno$'1'$data<- apply(HvDRR09$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'2'$data<- apply(HvDRR09$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'3'$data<- apply(HvDRR09$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'4'$data<- apply(HvDRR09$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'5'$data<- apply(HvDRR09$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'6'$data<- apply(HvDRR09$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'7'$data<- apply(HvDRR09$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR09)[1] <- "riself"

summary(HvDRR09)
plot(HvDRR09)
HvDRR09 <- subset(HvDRR09,ind=c(1:nrow(HvDRR09$pheno))[!is.na(HvDRR09$pheno)])

HvDRR09per <- calc.genoprob(HvDRR09, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR09out <- scanone(HvDRR09per, method = "hk")
plot(HvDRR09out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR09permout <- scanone(HvDRR09per, method = "hk", n.perm = 4000)
summary(HvDRR09out, perms = HvDRR09permout, alpha = 0.05, pvalues = TRUE)
HvDRR09qtl1 <- makeqtl(HvDRR09per, chr = 7, pos = 90, what = "prob")
HvDRR09outc5 <- addqtl(HvDRR09per, qtl = HvDRR09qtl1, method = "hk")
summary(HvDRR09outc5, perms = HvDRR09permout, alpha = 0.05, pvalues = TRUE)

HvDRR09rqtl <- refineqtl(HvDRR09per, qtl = HvDRR09qtl1, method = "hk", verbose = FALSE)
summary(HvDRR09rqtl)
summary(fitqtl(HvDRR09per, qtl = HvDRR09rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR09per, qtl = HvDRR09rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

find.marker(HvDRR09per, 7, 90)
effectplot(HvDRR09per, mname1 = "SCRI_RS_220780")
summary(HvDRR09per, mname1 = "SCRI_RS_220780")
plotPXG(HvDRR09per, "SCRI_RS_220780")

HvDRR09CI <- lodint(HvDRR09rqtl, chr = 7, qtl.index = 1)

pt(6.082, 98, lower.tail = FALSE)


#HvDRR10####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR10 <- Populations$HvDRR10

HvDRR10$geno$'1'$data<- apply(HvDRR10$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'2'$data<- apply(HvDRR10$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'3'$data<- apply(HvDRR10$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'4'$data<- apply(HvDRR10$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'5'$data<- apply(HvDRR10$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'6'$data<- apply(HvDRR10$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'7'$data<- apply(HvDRR10$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'1'$data<- apply(HvDRR10$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'2'$data<- apply(HvDRR10$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'3'$data<- apply(HvDRR10$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'4'$data<- apply(HvDRR10$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'5'$data<- apply(HvDRR10$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'6'$data<- apply(HvDRR10$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'7'$data<- apply(HvDRR10$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR10)[1] <- "riself"

summary(HvDRR10)
plot(HvDRR10)
HvDRR10 <- subset(HvDRR10,ind=c(1:nrow(HvDRR10$pheno))[!is.na(HvDRR10$pheno)])

HvDRR10per <- calc.genoprob(HvDRR10, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR10out <- scanone(HvDRR10per, method = "hk")
plot(HvDRR10out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR10permout <- scanone(HvDRR10per, method = "hk", n.perm = 10000)
summary(HvDRR10out, perms = HvDRR10permout, alpha = 0.05, pvalues = TRUE)
HvDRR10qtl1 <- makeqtl(HvDRR10per, chr = 3, pos = 124, what = "prob")
HvDRR10outc3 <- addqtl(HvDRR10per, qtl = HvDRR10qtl1, method = "hk")
summary(HvDRR10outc3, perms = HvDRR10permout, alpha = 0.05, pvalues = TRUE)

HvDRR10rqtl <- refineqtl(HvDRR10per, qtl = HvDRR10qtl1, method = "hk", verbose = FALSE)
summary(HvDRR10rqtl)
summary(fitqtl(HvDRR10per, qtl = HvDRR10rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR10per, qtl = HvDRR10rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

find.marker(HvDRR10per, 3, 123.7)
effectplot(HvDRR10per, mname1 = "JHI-Hv50k-2016-205404")
summary(HvDRR10per, mname1 = "JHI-Hv50k-2016-205404")
plotPXG(HvDRR10per, "JHI-Hv50k-2016-205404")

HvDRR10CI <- lodint(HvDRR10rqtl, chr = 3, qtl.index = 1, expandtomarkers = TRUE)

pt(10.03, 87, lower.tail = FALSE)
#p<0.01


#HvDRR11####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR11 <- Populations$HvDRR11

HvDRR11$geno$'1'$data<- apply(HvDRR11$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'2'$data<- apply(HvDRR11$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'3'$data<- apply(HvDRR11$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'4'$data<- apply(HvDRR11$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'5'$data<- apply(HvDRR11$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'6'$data<- apply(HvDRR11$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'7'$data<- apply(HvDRR11$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'1'$data<- apply(HvDRR11$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'2'$data<- apply(HvDRR11$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'3'$data<- apply(HvDRR11$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'4'$data<- apply(HvDRR11$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'5'$data<- apply(HvDRR11$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'6'$data<- apply(HvDRR11$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'7'$data<- apply(HvDRR11$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR11)[1] <- "riself"


summary(HvDRR11)
plot(HvDRR11)
HvDRR11 <- subset(HvDRR11,ind=c(1:nrow(HvDRR11$pheno))[!is.na(HvDRR11$pheno)])

testmarkername <- as.data.frame(HvDRR11$geno$`3`$map)
which(row.names(testmarkername) == 'SCRI_RS_150063')
marker2drop <- markernames(HvDRR11, chr = 3)[125]
HvDRR11 <- drop.markers(HvDRR11, marker2drop)

HvDRR11per <- calc.genoprob(HvDRR11, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR11out <- scanone(HvDRR11per, method = "hk")
plot(HvDRR11out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR11permout <- scanone(HvDRR11per, method = "hk", n.perm = 4000)
summary(HvDRR11out, perms = HvDRR11permout, alpha = 0.05, pvalues = TRUE)
HvDRR11qtl1 <- makeqtl(HvDRR11per, chr = c(2, 3), pos = c(166, 110), what = "prob")
HvDRR11outc23 <- addqtl(HvDRR11per, qtl = HvDRR11qtl1, method = "hk")
summary(HvDRR11outc23, perms = HvDRR11permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR11out, HvDRR11outc23, chr = c(-2, -3), col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR11out, HvDRR11outc23, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)
HvDRR11qtl2 <- makeqtl(HvDRR11per, c(2, 3, 7), c(166, 110, 53.4), what = "prob")
HvDRR11outc237 <- addqtl(HvDRR11per, qtl = HvDRR11qtl2, method = "hk")
summary(HvDRR11outc237, perms = HvDRR11permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR11out, HvDRR11outc237, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)

HvDRR11rqtl <- refineqtl(HvDRR11per, qtl = HvDRR11qtl2, method = "hk", verbose = FALSE)
summary(HvDRR11rqtl)
summary(fitqtl(HvDRR11per, qtl = HvDRR11rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR11per, 2, 166)
effectplot(HvDRR11per, mname1 = "JHI-Hv50k-2016-131359")
summary(HvDRR11per, mname1 = "JHI-Hv50k-2016-131359")
plotPXG(HvDRR11per, "JHI-Hv50k-2016-131359")

find.marker(HvDRR11per, 3, 110)
effectplot(HvDRR11per, mname1 = "JHI-Hv50k-2016-204079")
summary(HvDRR11per, mname1 = "JHI-Hv50k-2016-204079")
plotPXG(HvDRR11per, "JHI-Hv50k-2016-204079")

find.marker(HvDRR11per, 7, 53.4)
effectplot(HvDRR11per, mname1 = "JHI-Hv50k-2016-460580")
summary(HvDRR11per, mname1 = "JHI-Hv50k-2016-460580")
plotPXG(HvDRR11per, "JHI-Hv50k-2016-460580")

HvDRR11CI2 <- lodint(HvDRR11rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR11CI3 <- lodint(HvDRR11rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)
HvDRR11CI7 <- lodint(HvDRR11rqtl, chr = 7, qtl.index = 3, expandtomarkers = TRUE)

summary(fitqtl(HvDRR11per, qtl = HvDRR11rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3))

pt(-6.460, 93, lower.tail = FALSE)
#p<0.01
pt(6.950, 93, lower.tail = FALSE)
#p<0.01
pt(-4.366, 94, lower.tail = FALSE)
#p<0.01

#HvDRR12####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR12 <- Populations$HvDRR12

HvDRR12$geno$'1'$data<- apply(HvDRR12$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'2'$data<- apply(HvDRR12$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'3'$data<- apply(HvDRR12$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'4'$data<- apply(HvDRR12$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'5'$data<- apply(HvDRR12$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'6'$data<- apply(HvDRR12$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'7'$data<- apply(HvDRR12$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
#to eliminate heterozygous
HvDRR12$geno$'1'$data<- apply(HvDRR12$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'2'$data<- apply(HvDRR12$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'3'$data<- apply(HvDRR12$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'4'$data<- apply(HvDRR12$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'5'$data<- apply(HvDRR12$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'6'$data<- apply(HvDRR12$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'7'$data<- apply(HvDRR12$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR12)[1] <- "riself"

summary(HvDRR12)
plot(HvDRR12)
HvDRR12 <- subset(HvDRR12,ind=c(1:nrow(HvDRR12$pheno))[!is.na(HvDRR12$pheno)])


HvDRR12per <- calc.genoprob(HvDRR12, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR12out <- scanone(HvDRR12per, method = "hk")
plot(HvDRR12out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR12permout <- scanone(HvDRR12per, method = "hk", n.perm = 4000)
summary(HvDRR12out, perms = HvDRR12permout, alpha = 0.05, pvalues = TRUE)
HvDRR12qtl1 <- makeqtl(HvDRR12per, chr = c(5,7), pos = c(101, 41.5), what = "prob")
HvDRR12outc57 <- addqtl(HvDRR12per, qtl = HvDRR12qtl1, method = "hk")
summary(HvDRR12outc57, perms = HvDRR12permout, alpha = 0.05, pvalues = TRUE)

HvDRR12rqtl <- refineqtl(HvDRR12per, qtl = HvDRR12qtl1, method = "hk", verbose = FALSE)
summary(HvDRR12rqtl)
summary(fitqtl(HvDRR12per, qtl = HvDRR12rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR12per, qtl = HvDRR12rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1+Q2))

find.marker(HvDRR12per, 5, 8)
effectplot(HvDRR12per, mname1 = "JHI-Hv50k-2016-278416")
summary(HvDRR12per, mname1 = "JHI-Hv50k-2016-278416")
plotPXG(HvDRR12per, "JHI-Hv50k-2016-278416")

find.marker(HvDRR12per, 7, 41.5)
effectplot(HvDRR12per, mname1 = "BOPA2_12_30893")
summary(HvDRR12per, mname1 = "BOPA2_12_30893")
plotPXG(HvDRR12per, "BOPA2_12_30893")

HvDRR12CI5 <- lodint(HvDRR12rqtl, chr = 2, qtl.index = 1)
HvDRR12CI7 <- lodint(HvDRR12rqtl, chr = 7, qtl.index = 2)

pt(4.038,62, lower.tail = FALSE)
#p<0.01
pt(-8.375, 62, lower.tail = FALSE)
#p<0.01

#HvDRR13####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR13 <- Populations$HvDRR13

HvDRR13$geno$'1'$data<- apply(HvDRR13$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'2'$data<- apply(HvDRR13$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'3'$data<- apply(HvDRR13$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'4'$data<- apply(HvDRR13$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'5'$data<- apply(HvDRR13$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'6'$data<- apply(HvDRR13$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'7'$data<- apply(HvDRR13$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'1'$data<- apply(HvDRR13$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'2'$data<- apply(HvDRR13$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'3'$data<- apply(HvDRR13$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'4'$data<- apply(HvDRR13$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'5'$data<- apply(HvDRR13$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'6'$data<- apply(HvDRR13$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'7'$data<- apply(HvDRR13$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR13)[1] <- "riself"

summary(HvDRR13)
plot(HvDRR13)
HvDRR13 <- subset(HvDRR13,ind=c(1:nrow(HvDRR13$pheno))[!is.na(HvDRR13$pheno)])

HvDRR13per <- calc.genoprob(HvDRR13, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR13out <- scanone(HvDRR13per, method = "hk")
plot(HvDRR13out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR13permout <- scanone(HvDRR13per, method = "hk", n.perm = 4000)
summary(HvDRR13out, perms = HvDRR13permout, alpha = 0.05, pvalues = TRUE)
HvDRR13qtl1 <- makeqtl(HvDRR13per, chr = 7, pos = 55, what = "prob")
HvDRR13outc7 <- addqtl(HvDRR13per, qtl = HvDRR13qtl1, method = "hk")
summary(HvDRR13outc7, perms = HvDRR13permout, alpha = 0.05, pvalues = TRUE)

HvDRR13rqtl <- refineqtl(HvDRR13per, qtl = HvDRR13qtl1, method = "hk", verbose = FALSE)
summary(HvDRR13rqtl)
summary(fitqtl(HvDRR13per, qtl = HvDRR13rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR13per, qtl = HvDRR13rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

find.marker(HvDRR13per, 7, 55)
effectplot(HvDRR13per, mname1 = "BOPA2_12_10218")
summary(HvDRR13per, mname1 = "BOPA2_12_10218")
plotPXG(HvDRR13per, "BOPA2_12_10218")

HvDRR13CI7 <- lodint(HvDRR13rqtl, chr = 7, qtl.index = 1, expandtomarkers = TRUE)

pt(-11.7, 66, lower.tail = FALSE)
#p<0.01

#HvDRR14####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR14 <- Populations$HvDRR14

HvDRR14$geno$'1'$data<- apply(HvDRR14$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'2'$data<- apply(HvDRR14$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'3'$data<- apply(HvDRR14$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'4'$data<- apply(HvDRR14$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'5'$data<- apply(HvDRR14$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'6'$data<- apply(HvDRR14$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'7'$data<- apply(HvDRR14$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'1'$data<- apply(HvDRR14$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'2'$data<- apply(HvDRR14$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'3'$data<- apply(HvDRR14$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'4'$data<- apply(HvDRR14$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'5'$data<- apply(HvDRR14$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'6'$data<- apply(HvDRR14$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'7'$data<- apply(HvDRR14$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR14)[1] <- "riself"

summary(HvDRR14)
plot(HvDRR14)
HvDRR14 <- subset(HvDRR14,ind=c(1:nrow(HvDRR14$pheno))[!is.na(HvDRR14$pheno)])

HvDRR14per <- calc.genoprob(HvDRR14, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR14out <- scanone(HvDRR14per, method = "hk")
plot(HvDRR14out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR14permout <- scanone(HvDRR14per, method = "hk", n.perm = 4000)
summary(HvDRR14out, perms = HvDRR14permout, alpha = 0.05, pvalues = TRUE)
HvDRR14qtl1 <- makeqtl(HvDRR14per, chr = 7, pos = 71, what = "prob")
HvDRR14outc7 <- addqtl(HvDRR14per, qtl = HvDRR14qtl1, method = "hk")
summary(HvDRR14outc7, perms = HvDRR14permout, alpha = 0.05, pvalues = TRUE)

HvDRR14rqtl <- refineqtl(HvDRR14per, qtl = HvDRR14qtl1, method = "hk", verbose = FALSE)
summary(HvDRR14rqtl)
summary(fitqtl(HvDRR14per, qtl = HvDRR14rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR14per, 7, 71)
effectplot(HvDRR14per, mname1 = "SCRI_RS_121774")
summary(HvDRR14per, mname1 = "SCRI_RS_121774")
plotPXG(HvDRR14per, "SCRI_RS_121774")

summary(fitqtl(HvDRR14per, qtl = HvDRR14rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR14CI <- lodint(HvDRR14rqtl, chr = 7, qtl.index = 1, expandtomarkers = TRUE)

pt(-5.647, 74, lower.tail = FALSE)
#p<0.01

#HvDRR15####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR15 <- Populations$HvDRR15

HvDRR15$geno$'1'$data<- apply(HvDRR15$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'2'$data<- apply(HvDRR15$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'3'$data<- apply(HvDRR15$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'4'$data<- apply(HvDRR15$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'5'$data<- apply(HvDRR15$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'6'$data<- apply(HvDRR15$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'7'$data<- apply(HvDRR15$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'1'$data<- apply(HvDRR15$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'2'$data<- apply(HvDRR15$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'3'$data<- apply(HvDRR15$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'4'$data<- apply(HvDRR15$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'5'$data<- apply(HvDRR15$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'6'$data<- apply(HvDRR15$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'7'$data<- apply(HvDRR15$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR15)[1] <- "riself"

summary(HvDRR15)
plot(HvDRR15)
HvDRR15 <- subset(HvDRR15,ind=c(1:nrow(HvDRR15$pheno))[!is.na(HvDRR15$pheno)])

HvDRR15per <- calc.genoprob(HvDRR15, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR15out <- scanone(HvDRR15per, method = "hk")
plot(HvDRR15out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR15permout <- scanone(HvDRR15per, method = "hk", n.perm = 4000)
summary(HvDRR15out, perms = HvDRR15permout, alpha = 0.05, pvalues = TRUE)
HvDRR15qtl1 <- makeqtl(HvDRR15per, chr = 7, pos = 45, what = "prob")
HvDRR15outc7 <- addqtl(HvDRR15per, qtl = HvDRR15qtl1, method = "hk")
summary(HvDRR15outc7, perms = HvDRR15permout, alpha = 0.05, pvalues = TRUE)

HvDRR15rqtl <- refineqtl(HvDRR15per, qtl = HvDRR15qtl1, method = "hk", verbose = FALSE)
summary(HvDRR15rqtl)
summary(fitqtl(HvDRR15per, qtl = HvDRR15rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR15per, qtl = HvDRR15rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

find.marker(HvDRR15per, 7, 45)
effectplot(HvDRR15per, mname1 = "JHI-Hv50k-2016-460172")
summary(HvDRR15per, mname1 = "JHI-Hv50k-2016-460172")
plotPXG(HvDRR15per, "JHI-Hv50k-2016-460172")

HvDRR15CI <- lodint(HvDRR15rqtl, chr = 7, qtl.index = 1)

pt(-8.478, 45, lower.tail = FALSE)
#p<0.01

#HvDRR16####


Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR16 <- Populations$HvDRR16

HvDRR16$geno$'1'$data<- apply(HvDRR16$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'2'$data<- apply(HvDRR16$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'3'$data<- apply(HvDRR16$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'4'$data<- apply(HvDRR16$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'5'$data<- apply(HvDRR16$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'6'$data<- apply(HvDRR16$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'7'$data<- apply(HvDRR16$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'1'$data<- apply(HvDRR16$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'2'$data<- apply(HvDRR16$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'3'$data<- apply(HvDRR16$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'4'$data<- apply(HvDRR16$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'5'$data<- apply(HvDRR16$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'6'$data<- apply(HvDRR16$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'7'$data<- apply(HvDRR16$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR16)[1] <- "riself"

summary(HvDRR16)
plot(HvDRR16)
HvDRR16 <- subset(HvDRR16,ind=c(1:nrow(HvDRR16$pheno))[!is.na(HvDRR16$pheno)])

testmarkername <- as.data.frame(HvDRR16$geno$`2`$map)
which(row.names(testmarkername) == 'SCRI_RS_202469')
marker2drop <- markernames(HvDRR16, chr = 2)[318]
HvDRR16 <- drop.markers(HvDRR16, marker2drop)

HvDRR16per <- calc.genoprob(HvDRR16, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR16out <- scanone(HvDRR16per, method = "hk")
plot(HvDRR16out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR16permout <- scanone(HvDRR16per, method = "hk", n.perm = 4000)
summary(HvDRR16out, perms = HvDRR16permout, alpha = 0.05, pvalues = TRUE)
HvDRR16qtl1 <- makeqtl(HvDRR16per, chr = 2, pos = 50.2, what = "prob")
HvDRR16outc2 <- addqtl(HvDRR16per, qtl = HvDRR16qtl1, method = "hk")
summary(HvDRR16outc2, perms = HvDRR16permout, alpha = 0.05, pvalues = TRUE)
HvDRR16qtl2 <- makeqtl(HvDRR16per, c(2, 2), c(50.2,264), what = "prob")
HvDRR16outc22 <- addqtl(HvDRR16per, qtl = HvDRR16qtl2, method = "hk")
summary(HvDRR16outc22, perms = HvDRR16permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR16out, HvDRR16outc22, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)

HvDRR16rqtl <- refineqtl(HvDRR16per, qtl = HvDRR16qtl2, method = "hk", verbose = FALSE)
summary(HvDRR16rqtl)
summary(fitqtl(HvDRR16per, qtl = HvDRR16rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR16per, 2, 59.2)
effectplot(HvDRR16per, mname1 = "JHI-Hv50k-2016-73695")
summary(HvDRR16per, mname1 = "JHI-Hv50k-2016-73695")
plotPXG(HvDRR16per, "JHI-Hv50k-2016-73695")

find.marker(HvDRR16per, 2, 264.2)
effectplot(HvDRR16per, mname1 = "JHI-Hv50k-2016-129807")
summary(HvDRR16per, mname1 = "JHI-Hv50k-2016-129807")
plotPXG(HvDRR16per, "JHI-Hv50k-2016-129807")

summary(fitqtl(HvDRR16per, qtl = HvDRR16rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1+Q2))

HvDRR16CI2 <- lodint(HvDRR16rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR16CI22 <- lodint(HvDRR16rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)

pt(-12.37, 86, lower.tail = FALSE)
#p<0.01
pt(4.81, 86, lower.tail = FALSE)
#p<0.01




#HvDRR17####


Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR17 <- Populations$HvDRR17

HvDRR17$geno$'1'$data<- apply(HvDRR17$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'2'$data<- apply(HvDRR17$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'3'$data<- apply(HvDRR17$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'4'$data<- apply(HvDRR17$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'5'$data<- apply(HvDRR17$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'6'$data<- apply(HvDRR17$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'7'$data<- apply(HvDRR17$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'1'$data<- apply(HvDRR17$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'2'$data<- apply(HvDRR17$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'3'$data<- apply(HvDRR17$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'4'$data<- apply(HvDRR17$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'5'$data<- apply(HvDRR17$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'6'$data<- apply(HvDRR17$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'7'$data<- apply(HvDRR17$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR17)[1] <- "riself"

summary(HvDRR17)
plot(HvDRR17)
HvDRR17 <- subset(HvDRR17,ind=c(1:nrow(HvDRR17$pheno))[!is.na(HvDRR17$pheno)])

testmarkername <- as.data.frame(HvDRR17$geno$`1`$map)
which(row.names(testmarkername) == 'JHI-Hv50k-2016-55638')
marker2drop <- markernames(HvDRR17, chr = 1)[434]
HvDRR17 <- drop.markers(HvDRR17, marker2drop)

testmarkername <- as.data.frame(HvDRR17$geno$`2`$map)
which(row.names(testmarkername) == 'JHI-Hv50k-2016-73697')
marker2drop <- markernames(HvDRR17, chr = 2)[93]
HvDRR17 <- drop.markers(HvDRR17, marker2drop)

HvDRR17per <- calc.genoprob(HvDRR17, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR17out <- scanone(HvDRR17per, method = "hk")
plot(HvDRR17out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR17permout <- scanone(HvDRR17per, method = "hk", n.perm=4000)
summary(HvDRR17out, perms = HvDRR17permout, alpha = 0.05, pvalues = TRUE)
HvDRR17qtl1 <- makeqtl(HvDRR17per, chr = c(1,2), pos = c(49.3,53), what = "prob")
HvDRR17outc2 <- addqtl(HvDRR17per, qtl = HvDRR17qtl1, method = "hk")
summary(HvDRR17outc2, perms = HvDRR17permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR17out, HvDRR17outc2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
HvDRR17qtl2 <- makeqtl(HvDRR17per, c(1,2,2), c(49.3, 53, 108), what = "prob")
HvDRR17outc122 <- addqtl(HvDRR17per, qtl = HvDRR17qtl2, method = "hk")
summary(HvDRR17outc122, perms = HvDRR17permout, alpha = 0.05, pvalues = TRUE)

HvDRR17rqtl <- refineqtl(HvDRR17per, qtl = HvDRR17qtl2, method = "hk", verbose = FALSE)
summary(HvDRR17rqtl)
summary(fitqtl(HvDRR17per, qtl = HvDRR17rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR17per, 1, 222.6)
effectplot(HvDRR17per, mname1 = "JHI-Hv50k-2016-51526")
summary(HvDRR17per, mname1 = "JHI-Hv50k-2016-51526")
plotPXG(HvDRR17per, "JHI-Hv50k-2016-51526")

find.marker(HvDRR17per, 2, 53)
effectplot(HvDRR17per, mname1 = "JHI-Hv50k-2016-73780")
summary(HvDRR17per, mname1 = "JHI-Hv50k-2016-73780")
plotPXG(HvDRR17per, "JHI-Hv50k-2016-73780")

find.marker(HvDRR17per, 2, 107.5)
effectplot(HvDRR17per, mname1 = "JHI-Hv50k-2016-87930")
summary(HvDRR17per, mname1 = "JHI-Hv50k-2016-87930")
plotPXG(HvDRR17per, "JHI-Hv50k-2016-87930")

summary(fitqtl(HvDRR17per, qtl = HvDRR17rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1 + Q2 + Q3))

HvDRR17CI1 <- lodint(HvDRR17rqtl, chr = 1, qtl.index = 1, expandtomarkers = TRUE)
HvDRR17CI2 <- lodint(HvDRR17rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)
HvDRR17CI22 <- lodint(HvDRR17rqtl, chr = 2, qtl.index = 3, expandtomarkers = TRUE)

pt(2.572, 110, lower.tail = FALSE)
pt(-13.151, 110, lower.tail = FALSE)
pt(-5.353, 110, lower.tail = FALSE)


#HvDRR18####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR18 <- Populations$HvDRR18

HvDRR18$geno$'1'$data<- apply(HvDRR18$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'2'$data<- apply(HvDRR18$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'3'$data<- apply(HvDRR18$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'4'$data<- apply(HvDRR18$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'5'$data<- apply(HvDRR18$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'6'$data<- apply(HvDRR18$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'7'$data<- apply(HvDRR18$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'1'$data<- apply(HvDRR18$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'2'$data<- apply(HvDRR18$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'3'$data<- apply(HvDRR18$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'4'$data<- apply(HvDRR18$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'5'$data<- apply(HvDRR18$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'6'$data<- apply(HvDRR18$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'7'$data<- apply(HvDRR18$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR18)[1] <- "riself"

summary(HvDRR18)
plot(HvDRR18)
HvDRR18 <- subset(HvDRR18,ind=c(1:nrow(HvDRR18$pheno))[!is.na(HvDRR18$pheno)])

HvDRR18per <- calc.genoprob(HvDRR18, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR18out <- scanone(HvDRR18per, method = "hk")
plot(HvDRR18out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR18permout <- scanone(HvDRR18per, method = "hk", n.perm = 4000)
summary(HvDRR18out, perms = HvDRR18permout, alpha = 0.05, pvalues = TRUE)
HvDRR18qtl1 <- makeqtl(HvDRR18per, chr = 2, pos = 49, what = "prob")
HvDRR18outc2 <- addqtl(HvDRR18per, qtl = HvDRR18qtl1, method = "hk")
summary(HvDRR18outc2, perms = HvDRR18permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR18out, HvDRR18outc2, chr = -2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR18out, HvDRR18outc2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)

HvDRR18rqtl <- refineqtl(HvDRR18per, qtl = HvDRR18qtl1, method = "hk", verbose = FALSE)
summary(HvDRR18rqtl)
summary(fitqtl(HvDRR18per, qtl = HvDRR18rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR18per, 2, 49)
effectplot(HvDRR18per, mname1 = "JHI-Hv50k-2016-73692")
summary(HvDRR18per, mname1 = "JHI-Hv50k-2016-73692")
plotPXG(HvDRR18per, "JHI-Hv50k-2016-73692")

summary(fitqtl(HvDRR18per, qtl = HvDRR18rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1))

HvDRR18CI2 <- lodint(HvDRR18rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(7.286, 77, lower.tail = FALSE)
#p<0.01






























#HvDRR19####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR19 <- Populations$HvDRR19

HvDRR19$geno$'1'$data<- apply(HvDRR19$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'2'$data<- apply(HvDRR19$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'3'$data<- apply(HvDRR19$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'4'$data<- apply(HvDRR19$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'5'$data<- apply(HvDRR19$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'6'$data<- apply(HvDRR19$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'7'$data<- apply(HvDRR19$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'1'$data<- apply(HvDRR19$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'2'$data<- apply(HvDRR19$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'3'$data<- apply(HvDRR19$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'4'$data<- apply(HvDRR19$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'5'$data<- apply(HvDRR19$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'6'$data<- apply(HvDRR19$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'7'$data<- apply(HvDRR19$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR19)[1] <- "riself"

summary(HvDRR19)
plot(HvDRR19)
HvDRR19 <- subset(HvDRR19,ind=c(1:nrow(HvDRR19$pheno))[!is.na(HvDRR19$pheno)])

HvDRR19per <- calc.genoprob(HvDRR19, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR19out <- scanone(HvDRR19per, method = "hk")
plot(HvDRR19out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR19permout <- scanone(HvDRR19per, method = "hk", n.perm = 4000)
summary(HvDRR19out, perms = HvDRR19permout, alpha = 0.05, pvalues = TRUE)
HvDRR19qtl1 <- makeqtl(HvDRR19per, chr = 2, pos = 100, what = "prob")
HvDRR19outc2 <- addqtl(HvDRR19per, qtl = HvDRR19qtl1, method = "hk")
summary(HvDRR19outc2, perms = HvDRR19permout, alpha = 0.05, pvalues = TRUE)

HvDRR19rqtl <- refineqtl(HvDRR19per, qtl = HvDRR19qtl1, method = "hk", verbose = FALSE)
summary(HvDRR19rqtl)
summary(fitqtl(HvDRR19per, qtl = HvDRR19rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR19per, 2, 100)
effectplot(HvDRR19per, mname1 = "BOPA1_ABC08774-1-1-752")
summary(HvDRR19per, mname1 = "BOPA1_ABC08774-1-1-752")
plotPXG(HvDRR19per, "BOPA1_ABC08774-1-1-752")

summary(fitqtl(HvDRR19per, qtl = HvDRR19rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1))

HvDRR19CI2 <- lodint(HvDRR19rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(-6.311, 81, lower.tail = FALSE)
#p<0.01

















#HvDRR20####

#setwd("/gpfs/project/projects/qggp/Francesco_barley/1stpaper/Rprojects/QTLanalysis/")
Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
#library(qtl)
HvDRR20 <- Populations$HvDRR20

HvDRR20$geno$'1'$data<- apply(HvDRR20$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'2'$data<- apply(HvDRR20$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'3'$data<- apply(HvDRR20$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'4'$data<- apply(HvDRR20$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'5'$data<- apply(HvDRR20$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'6'$data<- apply(HvDRR20$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'7'$data<- apply(HvDRR20$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
#to eliminate heterozygous
HvDRR20$geno$'1'$data<- apply(HvDRR20$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'2'$data<- apply(HvDRR20$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'3'$data<- apply(HvDRR20$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'4'$data<- apply(HvDRR20$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'5'$data<- apply(HvDRR20$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'6'$data<- apply(HvDRR20$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'7'$data<- apply(HvDRR20$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR20)[1] <- "riself"

summary(HvDRR20)
plot(HvDRR20)
HvDRR20 <- subset(HvDRR20,ind=c(1:nrow(HvDRR20$pheno))[!is.na(HvDRR20$pheno)])

HvDRR20per <- calc.genoprob(HvDRR20, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR20out <- scanone(HvDRR20per, method = "hk")
plot(HvDRR20out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR20permout <- scanone(HvDRR20per, method = "hk", n.perm = 4000)
summary(HvDRR20out, perms = HvDRR20permout, alpha = 0.05, pvalues = TRUE)
HvDRR20qtl1 <- makeqtl(HvDRR20per, chr = 2, pos = 117, what = "prob")
HvDRR20outc2 <- addqtl(HvDRR20per, qtl = HvDRR20qtl1, method = "hk")
summary(HvDRR20outc2, perms = HvDRR20permout, alpha = 0.05, pvalues = TRUE)
HvDRR20qtl2 <- makeqtl(HvDRR20per, chr = c(2,2), pos = c(117,30), what = "prob")
HvDRR20outc22 <- addqtl(HvDRR20per, qtl = HvDRR20qtl2, method = "hk")
summary(HvDRR20outc22, perms = HvDRR20permout, alpha = 0.05, pvalues = TRUE)
HvDRR20qtl3 <- makeqtl(HvDRR20per, chr = c(2,2,7), pos = c(117,30,43), what = "prob")
HvDRR20outc227 <- addqtl(HvDRR20per, qtl = HvDRR20qtl3, method = "hk")
summary(HvDRR20outc227, perms = HvDRR20permout, alpha = 0.05, pvalues = TRUE)

HvDRR20rqtl <- refineqtl(HvDRR20per, qtl = HvDRR20qtl3, method = "hk", verbose = FALSE)
summary(HvDRR20rqtl)
summary(fitqtl(HvDRR20per, qtl = HvDRR20rqtl, method = "hk"), pvalues = FALSE)


find.marker(HvDRR20per, 2, 114.9)
effectplot(HvDRR20per, mname1 = "BOPA1_ABC08774-1-1-752")
summary(HvDRR20per, mname1 = "BOPA1_ABC08774-1-1-752")
plotPXG(HvDRR20per, "BOPA1_ABC08774-1-1-752")

find.marker(HvDRR20per, 2, 33)
effectplot(HvDRR20per, mname1 = "JHI-Hv50k-2016-73370")
summary(HvDRR20per, mname1 = "JHI-Hv50k-2016-73370")
plotPXG(HvDRR20per, "JHI-Hv50k-2016-73370")

find.marker(HvDRR20per, 7, 43)
effectplot(HvDRR20per, mname1 = "BOPA2_12_30894")
summary(HvDRR20per, mname1 = "BOPA2_12_30894")
plotPXG(HvDRR20per, "BOPA2_12_30894")

summary(fitqtl(HvDRR20per, qtl = HvDRR20rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3))

HvDRR20CI2 <- lodint(HvDRR20rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR20CI22 <- lodint(HvDRR20rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)
HvDRR20CI7 <- lodint(HvDRR20rqtl, chr = 7, qtl.index = 3, expandtomarkers = TRUE)

pt(6.512, 94, lower.tail = FALSE)
#p<0.01
pt(-6.037, 94, lower.tail = FALSE)
#p<0.01
pt(4.421, 94, lower.tail = FALSE)
#p<0.01













#HvDRR21####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR21 <- Populations$HvDRR21

HvDRR21$geno$'1'$data<- apply(HvDRR21$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'2'$data<- apply(HvDRR21$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'3'$data<- apply(HvDRR21$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'4'$data<- apply(HvDRR21$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'5'$data<- apply(HvDRR21$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'6'$data<- apply(HvDRR21$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'7'$data<- apply(HvDRR21$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'1'$data<- apply(HvDRR21$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'2'$data<- apply(HvDRR21$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'3'$data<- apply(HvDRR21$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'4'$data<- apply(HvDRR21$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'5'$data<- apply(HvDRR21$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'6'$data<- apply(HvDRR21$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'7'$data<- apply(HvDRR21$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR21)[1] <- "riself"

summary(HvDRR21)
plot(HvDRR21)
HvDRR21 <- subset(HvDRR21,ind=c(1:nrow(HvDRR21$pheno))[!is.na(HvDRR21$pheno)])
HvDRR21per <- calc.genoprob(HvDRR21, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR21out <- scanone(HvDRR21per, method = "hk")
plot(HvDRR21out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR21permout <- scanone(HvDRR21per, method = "hk", n.perm = 4000)
summary(HvDRR21out, perms = HvDRR21permout, alpha = 0.05, pvalues = TRUE)


#HvDRR22####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR22 <- Populations$HvDRR22

HvDRR22$geno$'1'$data<- apply(HvDRR22$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'2'$data<- apply(HvDRR22$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'3'$data<- apply(HvDRR22$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'4'$data<- apply(HvDRR22$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'5'$data<- apply(HvDRR22$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'6'$data<- apply(HvDRR22$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'7'$data<- apply(HvDRR22$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'1'$data<- apply(HvDRR22$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'2'$data<- apply(HvDRR22$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'3'$data<- apply(HvDRR22$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'4'$data<- apply(HvDRR22$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'5'$data<- apply(HvDRR22$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'6'$data<- apply(HvDRR22$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'7'$data<- apply(HvDRR22$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR22)[1] <- "riself"

summary(HvDRR22)
plot(HvDRR22)
HvDRR22 <- subset(HvDRR22,ind=c(1:nrow(HvDRR22$pheno))[!is.na(HvDRR22$pheno)])

HvDRR22per <- calc.genoprob(HvDRR22, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR22out <- scanone(HvDRR22per, method = "hk")
plot(HvDRR22out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR22permout <- scanone(HvDRR22per, method = "hk", n.perm=4000)
summary(HvDRR22out, perms = HvDRR22permout, alpha = 0.05, pvalues = TRUE)

HvDRR22qtl1 <- makeqtl(HvDRR22per, chr = c(2,5,7), pos = c(84.7,77.8,62), what = "prob")
HvDRR22outc7 <- addqtl(HvDRR22per, qtl = HvDRR22qtl1, method = "hk")
summary(HvDRR22outc7, perms = HvDRR22permout, alpha = 0.05, pvalues = TRUE)

HvDRR22rqtl <- refineqtl(HvDRR22per, qtl = HvDRR22qtl1, method = "hk", verbose = FALSE)
summary(HvDRR22rqtl)
summary(fitqtl(HvDRR22per, qtl = HvDRR22rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR22per, 2, 88.1)
effectplot(HvDRR22per, mname1 = "JHI-Hv50k-2016-90437")
summary(HvDRR22per, mname1 = "JHI-Hv50k-2016-90437")
plotPXG(HvDRR22per, "JHI-Hv50k-2016-90437")

find.marker(HvDRR22per, 5, 79.2)
effectplot(HvDRR22per, mname1 = "SCRI_RS_223712")
summary(HvDRR22per, mname1 = "SCRI_RS_223712")
plotPXG(HvDRR22per, "SCRI_RS_223712")

find.marker(HvDRR22per, 7, 63)
effectplot(HvDRR22per, mname1 = "JHI-Hv50k-2016-460104")
summary(HvDRR22per, mname1 = "JHI-Hv50k-2016-460104")
plotPXG(HvDRR22per, "JHI-Hv50k-2016-460104")

summary(fitqtl(HvDRR22per, qtl = HvDRR22rqtl, method = "hk", get.ests = TRUE, formula = Flowering ~Q1 + Q2 + Q3))

HvDRR22CI2 <- lodint(HvDRR22rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR22CI5 <- lodint(HvDRR22rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)
HvDRR22CI7 <- lodint(HvDRR22rqtl, chr = 7, qtl.index = 3, expandtomarkers = TRUE)

pt(-6.418, 86, lower.tail = FALSE)
#p<0.01
pt(-4.459, 86, lower.tail = FALSE)
#p<0.01
pt(-10.572, 86, lower.tail = FALSE)
#p<0.01



#HvDRR23####


Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR23 <- Populations$HvDRR23

HvDRR23$geno$'1'$data<- apply(HvDRR23$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'2'$data<- apply(HvDRR23$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'3'$data<- apply(HvDRR23$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'4'$data<- apply(HvDRR23$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'5'$data<- apply(HvDRR23$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'6'$data<- apply(HvDRR23$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'7'$data<- apply(HvDRR23$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'1'$data<- apply(HvDRR23$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'2'$data<- apply(HvDRR23$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'3'$data<- apply(HvDRR23$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'4'$data<- apply(HvDRR23$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'5'$data<- apply(HvDRR23$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'6'$data<- apply(HvDRR23$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'7'$data<- apply(HvDRR23$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR23)[1] <- "riself"

summary(HvDRR23)
plot(HvDRR23)
HvDRR23 <- subset(HvDRR23,ind=c(1:nrow(HvDRR23$pheno))[!is.na(HvDRR23$pheno)])

HvDRR23per <- calc.genoprob(HvDRR23, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR23out <- scanone(HvDRR23per, method = "hk")
plot(HvDRR23out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR23permout <- scanone(HvDRR23per, method = "hk", n.perm = 4000)
summary(HvDRR23out, perms = HvDRR23permout, alpha = 0.05, pvalues = TRUE)
HvDRR23qtl1 <- makeqtl(HvDRR23per, chr = 2, pos = 39, what = "prob")
HvDRR23outc2 <- addqtl(HvDRR23per, qtl = HvDRR23qtl1, method = "hk")
summary(HvDRR23outc2, perms = HvDRR23permout, alpha = 0.05, pvalues = TRUE)

HvDRR23rqtl <- refineqtl(HvDRR23per, qtl = HvDRR23qtl1, method = "hk", verbose = FALSE)
summary(HvDRR23rqtl)
summary(fitqtl(HvDRR23per, qtl = HvDRR23rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR23per, 2, 39)
effectplot(HvDRR23per, mname1 = "BOPA2_12_30870")
summary(HvDRR23per, mname1 = "BOPA2_12_30870")
plotPXG(HvDRR23per, "BOPA2_12_30870")

summary(fitqtl(HvDRR23per, qtl = HvDRR23rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1))

HvDRR23CI2 <- lodint(HvDRR23rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(-4.135, 84, lower.tail = FALSE)
#p<0.01


#HvDRR24####


Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR24 <- Populations$HvDRR24

HvDRR24$geno$'1'$data<- apply(HvDRR24$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'2'$data<- apply(HvDRR24$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'3'$data<- apply(HvDRR24$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'4'$data<- apply(HvDRR24$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'5'$data<- apply(HvDRR24$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'6'$data<- apply(HvDRR24$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'7'$data<- apply(HvDRR24$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'1'$data<- apply(HvDRR24$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'2'$data<- apply(HvDRR24$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'3'$data<- apply(HvDRR24$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'4'$data<- apply(HvDRR24$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'5'$data<- apply(HvDRR24$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'6'$data<- apply(HvDRR24$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'7'$data<- apply(HvDRR24$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR24)[1] <- "riself"

summary(HvDRR24)
plot(HvDRR24)
HvDRR24 <- subset(HvDRR24,ind=c(1:nrow(HvDRR24$pheno))[!is.na(HvDRR24$pheno)])

HvDRR24per <- calc.genoprob(HvDRR24, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR24out <- scanone(HvDRR24per, method = "hk")
plot(HvDRR24out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR24permout <- scanone(HvDRR24per, method = "hk", n.perm=4000)
summary(HvDRR24out, perms = HvDRR24permout, alpha = 0.05, pvalues = TRUE)
HvDRR24qtl1 <- makeqtl(HvDRR24per, chr = 2, pos = 97, what = "prob")
HvDRR24outc2 <- addqtl(HvDRR24per, qtl = HvDRR24qtl1, method = "hk")
summary(HvDRR24outc2, perms = HvDRR24permout, alpha = 0.05, pvalues = TRUE)

HvDRR24rqtl <- refineqtl(HvDRR24per, qtl = HvDRR24qtl1, method = "hk", verbose = FALSE)
summary(HvDRR24rqtl)
summary(fitqtl(HvDRR24per, qtl = HvDRR24rqtl, method = "hk"), pvalues = FALSE)


find.marker(HvDRR24per, 2, 97)
effectplot(HvDRR24per, mname1 = "JHI-Hv50k-2016-87488")
summary(HvDRR24per, mname1 = "JHI-Hv50k-2016-87488")
plotPXG(HvDRR24per, "JHI-Hv50k-2016-87488")

summary(fitqtl(HvDRR24per, qtl = HvDRR24rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR24CI2 <- lodint(HvDRR24rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(4.384, 62, lower.tail = FALSE)
#p<0.01









#HvDRR25####


Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR25 <- Populations$HvDRR25

HvDRR25$geno$'1'$data<- apply(HvDRR25$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'2'$data<- apply(HvDRR25$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'3'$data<- apply(HvDRR25$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'4'$data<- apply(HvDRR25$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'5'$data<- apply(HvDRR25$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'6'$data<- apply(HvDRR25$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'7'$data<- apply(HvDRR25$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'1'$data<- apply(HvDRR25$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'2'$data<- apply(HvDRR25$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'3'$data<- apply(HvDRR25$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'4'$data<- apply(HvDRR25$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'5'$data<- apply(HvDRR25$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'6'$data<- apply(HvDRR25$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'7'$data<- apply(HvDRR25$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR25)[1] <- "riself"
summary(HvDRR25)
plot(HvDRR25)
HvDRR25 <- subset(HvDRR25,ind=c(1:nrow(HvDRR25$pheno))[!is.na(HvDRR25$pheno)])

HvDRR25per <- calc.genoprob(HvDRR25, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR25out <- scanone(HvDRR25per, method = "hk")
plot(HvDRR25out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR25permout <- scanone(HvDRR25per, method = "hk", n.perm=10000)
summary(HvDRR25out, perms = HvDRR25permout, alpha = 0.05, pvalues = TRUE)
HvDRR25qtl1 <- makeqtl(HvDRR25per, chr = 2, pos = 45, what = "prob")
HvDRR25outc2 <- addqtl(HvDRR25per, qtl = HvDRR25qtl1, method = "hk")
summary(HvDRR25outc2, perms = HvDRR25permout, alpha = 0.05, pvalues = TRUE)
HvDRR25qtl2 <- makeqtl(HvDRR25per, c(2,6), c(45,75.6), what = "prob")
HvDRR25outc26 <- addqtl(HvDRR25per, qtl = HvDRR25qtl2, method = "hk")
summary(HvDRR25outc26, perms = HvDRR25permout, alpha = 0.05, pvalues = TRUE)

HvDRR25rqtl <- refineqtl(HvDRR25per, qtl = HvDRR25qtl2, method = "hk", verbose = FALSE)
summary(HvDRR25rqtl)
summary(fitqtl(HvDRR25per, qtl = HvDRR25rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR25per, 2,45)
effectplot(HvDRR25per, mname1 = "BK_14")
summary(HvDRR25per, mname1 = "BK_14")
plotPXG(HvDRR25per, "BK_14")

find.marker(HvDRR25per, 6,75.6)
effectplot(HvDRR25per, mname1 = "JHI-Hv50k-2016-382481")
summary(HvDRR25per, mname1 = "JHI-Hv50k-2016-382481")
plotPXG(HvDRR25per, "JHI-Hv50k-2016-382481")
summary(fitqtl(HvDRR25per, qtl = HvDRR25rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2))

HvDRR25CI2 <- lodint(HvDRR25rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR25CI6 <- lodint(HvDRR25rqtl, chr = 6, qtl.index = 2, expandtomarkers = TRUE)

pt(-11.221, 72, lower.tail = FALSE)
#p<0.01
pt(4.216, 72, lower.tail = FALSE)
#p<0.01





#HvDRR26####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR26 <- Populations$HvDRR26

HvDRR26$geno$'1'$data<- apply(HvDRR26$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'2'$data<- apply(HvDRR26$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'3'$data<- apply(HvDRR26$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'4'$data<- apply(HvDRR26$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'5'$data<- apply(HvDRR26$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'6'$data<- apply(HvDRR26$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'7'$data<- apply(HvDRR26$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'1'$data<- apply(HvDRR26$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'2'$data<- apply(HvDRR26$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'3'$data<- apply(HvDRR26$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'4'$data<- apply(HvDRR26$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'5'$data<- apply(HvDRR26$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'6'$data<- apply(HvDRR26$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'7'$data<- apply(HvDRR26$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR26)[1] <- "riself"

summary(HvDRR26)
plot(HvDRR26)
HvDRR26 <- subset(HvDRR26,ind=c(1:nrow(HvDRR26$pheno))[!is.na(HvDRR26$pheno)])


HvDRR26per <- calc.genoprob(HvDRR26, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR26out <- scanone(HvDRR26per, method = "hk")
plot(HvDRR26out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR26permout <- scanone(HvDRR26per, method = "hk", n.perm=4000)
summary(HvDRR26out, perms = HvDRR26permout, alpha = 0.05, pvalues = TRUE)

HvDRR26qtl1 <- makeqtl(HvDRR26per, chr = c(2,7), pos = c(35.1, 156), what = "prob")
HvDRR26outc27 <- addqtl(HvDRR26per, qtl = HvDRR26qtl1, method = "hk")
summary(HvDRR26outc27, perms = HvDRR26permout, alpha = 0.05, pvalues = TRUE)
HvDRR26qtl2 <- makeqtl(HvDRR26per, c(2,7,2), c(35.1, 156, 70), what = "prob")
HvDRR26outc272 <- addqtl(HvDRR26per, qtl = HvDRR26qtl2, method = "hk")
summary(HvDRR26outc272, perms = HvDRR26permout, alpha = 0.05, pvalues = TRUE)

HvDRR26rqtl <- refineqtl(HvDRR26per, qtl = HvDRR26qtl2, method = "hk", verbose = FALSE)
summary(HvDRR26rqtl)
summary(fitqtl(HvDRR26per, qtl = HvDRR26rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR26per, 2, 34)
effectplot(HvDRR26per, mname1 = "JHI-Hv50k-2016-73562")
summary(HvDRR26per, mname1 = "JHI-Hv50k-2016-73562")
plotPXG(HvDRR26per, "JHI-Hv50k-2016-73562")

find.marker(HvDRR26per,7, 157)
effectplot(HvDRR26per, mname1 = "JHI-Hv50k-2016-497524")
summary(HvDRR26per, mname1 = "JHI-Hv50k-2016-497524")
plotPXG(HvDRR26per, "JHI-Hv50k-2016-497524")

find.marker(HvDRR26per, 2, 70)
effectplot(HvDRR26per, mname1 = "BOPA1_6846-1161")
summary(HvDRR26per, mname1 = "BOPA1_6846-1161")
plotPXG(HvDRR26per, "BOPA1_6846-1161")

summary(fitqtl(HvDRR26per, qtl = HvDRR26rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3))

HvDRR26CI2 <- lodint(HvDRR26rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR26CI7 <- lodint(HvDRR26rqtl, chr = 7, qtl.index = 2, expandtomarkers = TRUE)
HvDRR26CI22 <- lodint(HvDRR26rqtl, chr = 2, qtl.index = 3, expandtomarkers = TRUE)

pt(4.763, 50, lower.tail = FALSE)
#p<0.01
pt(-3.417, 50, lower.tail = FALSE)
#p<0.01
pt(4.366, 50, lower.tail = FALSE)
#p<0.01



#HvDRR27####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR27 <- Populations$HvDRR27

HvDRR27$geno$'1'$data<- apply(HvDRR27$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'2'$data<- apply(HvDRR27$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'3'$data<- apply(HvDRR27$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'4'$data<- apply(HvDRR27$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'5'$data<- apply(HvDRR27$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'6'$data<- apply(HvDRR27$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'7'$data<- apply(HvDRR27$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'1'$data<- apply(HvDRR27$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'2'$data<- apply(HvDRR27$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'3'$data<- apply(HvDRR27$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'4'$data<- apply(HvDRR27$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'5'$data<- apply(HvDRR27$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'6'$data<- apply(HvDRR27$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'7'$data<- apply(HvDRR27$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR27)[1] <- "riself"

summary(HvDRR27)
plot(HvDRR27)
HvDRR27 <- subset(HvDRR27,ind=c(1:nrow(HvDRR27$pheno))[!is.na(HvDRR27$pheno)])

HvDRR27per <- calc.genoprob(HvDRR27, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR27out <- scanone(HvDRR27per, method = "hk")
plot(HvDRR27out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR27permout <- scanone(HvDRR27per, method = "hk", n.perm=4000)
summary(HvDRR27out, perms = HvDRR27permout, alpha = 0.05, pvalues = TRUE)
HvDRR27qtl1 <- makeqtl(HvDRR27per, chr = 2, pos = 44.4, what = "prob")
HvDRR27outc2 <- addqtl(HvDRR27per, qtl = HvDRR27qtl1, method = "hk")
summary(HvDRR27outc2, perms = HvDRR27permout, alpha = 0.05, pvalues = TRUE)

HvDRR27rqtl <- refineqtl(HvDRR27per, qtl = HvDRR27qtl1, method = "hk", verbose = FALSE)
summary(HvDRR27rqtl)
summary(fitqtl(HvDRR27per, qtl = HvDRR27rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR27per, qtl = HvDRR27rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

find.marker(HvDRR27per, 2, 44.4)
effectplot(HvDRR27per, mname1 = "BOPA2_12_30870")
summary(HvDRR27per, mname1 = "BOPA2_12_30870")
plotPXG(HvDRR27per, "BOPA2_12_30870")

HvDRR27CI2 <- lodint(HvDRR27rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(17.83, 91, lower.tail = FALSE)
#p<0.01


#HvDRR28####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR28 <- Populations$HvDRR28

HvDRR28$geno$'1'$data<- apply(HvDRR28$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'2'$data<- apply(HvDRR28$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'3'$data<- apply(HvDRR28$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'4'$data<- apply(HvDRR28$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'5'$data<- apply(HvDRR28$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'6'$data<- apply(HvDRR28$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'7'$data<- apply(HvDRR28$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'1'$data<- apply(HvDRR28$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'2'$data<- apply(HvDRR28$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'3'$data<- apply(HvDRR28$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'4'$data<- apply(HvDRR28$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'5'$data<- apply(HvDRR28$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'6'$data<- apply(HvDRR28$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'7'$data<- apply(HvDRR28$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR28)[1] <- "riself"

summary(HvDRR28)
plot(HvDRR28)
HvDRR28 <- subset(HvDRR28,ind=c(1:nrow(HvDRR28$pheno))[!is.na(HvDRR28$pheno)])

HvDRR28per <- calc.genoprob(HvDRR28, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR28out <- scanone(HvDRR28per, method = "hk")
plot(HvDRR28out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR28permout <- scanone(HvDRR28per, method = "hk", n.perm=4000)
summary(HvDRR28out, perms = HvDRR28permout, alpha = 0.05, pvalues = TRUE)
HvDRR28qtl1 <- makeqtl(HvDRR28per, chr = 2, pos = 48.3, what = "prob")
HvDRR28outc2 <- addqtl(HvDRR28per, qtl = HvDRR28qtl1, method = "hk")
summary(HvDRR28outc2, perms = HvDRR28permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR28out, HvDRR28outc2, chr = -2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR28out, HvDRR28outc2, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)
HvDRR28qtl2 <- makeqtl(HvDRR28per, c(2,2), c(48.3, 123), what = "prob")
HvDRR28outc22 <- addqtl(HvDRR28per, qtl = HvDRR28qtl2, method = "hk")
summary(HvDRR28outc22, perms = HvDRR28permout, alpha = 0.05, pvalues = TRUE)
HvDRR28qtl3 <- makeqtl(HvDRR28per, c(2,2,5), c(48.3, 123, 186), what = "prob")
HvDRR28outc225 <- addqtl(HvDRR28per, qtl = HvDRR28qtl3, method = "hk")
summary(HvDRR28outc225, perms = HvDRR28permout, alpha = 0.05, pvalues = TRUE)
HvDRR28qtl4 <- makeqtl(HvDRR28per, c(2,2,5,4), c(48.3, 123, 186, 70), what = "prob")
HvDRR28outc2257 <- addqtl(HvDRR28per, qtl = HvDRR28qtl4, method = "hk")
summary(HvDRR28outc2257, perms = HvDRR28permout, alpha = 0.05, pvalues = TRUE)

HvDRR28rqtl <- refineqtl(HvDRR28per, qtl = HvDRR28qtl4, method = "hk", verbose = FALSE)
summary(HvDRR28rqtl)
summary(fitqtl(HvDRR28per, qtl = HvDRR28rqtl, method = "hk"), pvalues = FALSE)


find.marker(HvDRR28per, 2, 48)
effectplot(HvDRR28per, mname1 = "JHI-Hv50k-2016-74702")
summary(HvDRR28per, mname1 = "JHI-Hv50k-2016-74702")
plotPXG(HvDRR28per, "JHI-Hv50k-2016-74702")

find.marker(HvDRR28per, 2, 123.2)
effectplot(HvDRR28per, mname1 = "BOPA1_4136-869")
summary(HvDRR28per, mname1 = "BOPA1_4136-869")
plotPXG(HvDRR28per, "BOPA1_4136-869")

find.marker(HvDRR28per, 5, 187)
effectplot(HvDRR28per, mname1 = "JHI-Hv50k-2016-335896")
summary(HvDRR28per, mname1 = "JHI-Hv50k-2016-335896")
plotPXG(HvDRR28per, "JHI-Hv50k-2016-335896")

find.marker(HvDRR28per, 4,70)
effectplot(HvDRR28per, mname1 = "JHI-Hv50k-2016-250940")
summary(HvDRR28per, mname1 = "JHI-Hv50k-2016-250940")
plotPXG(HvDRR28per, "JHI-Hv50k-2016-250940")

summary(fitqtl(HvDRR28per, qtl = HvDRR28rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1+Q2+Q3+Q4))

HvDRR28CI2 <- lodint(HvDRR28rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR28CI22 <- lodint(HvDRR28rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)
HvDRR28CI5 <- lodint(HvDRR28rqtl, chr = 5, qtl.index = 3, expandtomarkers = TRUE)
HvDRR28CI4 <- lodint(HvDRR28rqtl, chr = 4, qtl.index = 4, expandtomarkers = TRUE)

pt(5.899, 76, lower.tail = FALSE)
#p<0.01
pt(6.371, 76, lower.tail = FALSE)
#p<0.01
pt(-4.308, 76, lower.tail = FALSE)
#p<0.01
pt(-4.106, 76, lower.tail = FALSE)
#p<0.01



#HvDRR29####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR29 <- Populations$HvDRR29

HvDRR29$geno$'1'$data<- apply(HvDRR29$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'2'$data<- apply(HvDRR29$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'3'$data<- apply(HvDRR29$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'4'$data<- apply(HvDRR29$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'5'$data<- apply(HvDRR29$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'6'$data<- apply(HvDRR29$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'7'$data<- apply(HvDRR29$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'1'$data<- apply(HvDRR29$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'2'$data<- apply(HvDRR29$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'3'$data<- apply(HvDRR29$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'4'$data<- apply(HvDRR29$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'5'$data<- apply(HvDRR29$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'6'$data<- apply(HvDRR29$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'7'$data<- apply(HvDRR29$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
class(HvDRR29)[1] <- "riself"

summary(HvDRR29)
plot(HvDRR29)
HvDRR29 <- subset(HvDRR29,ind=c(1:nrow(HvDRR29$pheno))[!is.na(HvDRR29$pheno)])

HvDRR29per <- calc.genoprob(HvDRR29, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR29out <- scanone(HvDRR29per, method = "hk")
plot(HvDRR29out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR29permout <- scanone(HvDRR29per, method = "hk", n.perm=4000)
summary(HvDRR29out, perms = HvDRR29permout, alpha = 0.05, pvalues = TRUE)
HvDRR29qtl1 <- makeqtl(HvDRR29per, chr = 2, pos = 43, what = "prob")
HvDRR29outc2 <- addqtl(HvDRR29per, qtl = HvDRR29qtl1, method = "hk")
summary(HvDRR29outc2, perms = HvDRR29permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR29out, HvDRR29outc2, chr = -2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR29out, HvDRR29outc2, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)

HvDRR29qtl2 <- makeqtl(HvDRR29per, c(2,2,4), c(43, 110, 170), what = "prob")
HvDRR29outc224 <- addqtl(HvDRR29per, qtl = HvDRR29qtl2, method = "hk")
summary(HvDRR29outc224, perms = HvDRR29permout, alpha = 0.05, pvalues = TRUE)
HvDRR29qtl3 <- makeqtl(HvDRR29per, c(2,2,4,2), c(43, 110, 170, 207), what = "prob")
HvDRR29outc2242 <- addqtl(HvDRR29per, qtl = HvDRR29qtl3, method = "hk")
summary(HvDRR29outc2242, perms = HvDRR29permout, alpha = 0.05, pvalues = TRUE)

HvDRR29rqtl <- refineqtl(HvDRR29per, qtl = HvDRR29qtl3, method = "hk", verbose = FALSE)
summary(HvDRR29rqtl)
summary(fitqtl(HvDRR29per, qtl = HvDRR29rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR29per, 2, 42)
effectplot(HvDRR29per, mname1 = "JHI-Hv50k-2016-73584")
summary(HvDRR29per, mname1 = "JHI-Hv50k-2016-73584")
plotPXG(HvDRR29per, "JHI-Hv50k-2016-73584")

find.marker(HvDRR29per, 2, 109)
effectplot(HvDRR29per, mname1 = "JHI-Hv50k-2016-94652")
summary(HvDRR29per, mname1 = "JHI-Hv50k-2016-94652")
plotPXG(HvDRR29per, "JHI-Hv50k-2016-94652")

find.marker(HvDRR29per, 4,169.8)
effectplot(HvDRR29per, mname1 = "JHI-Hv50k-2016-272251")
summary(HvDRR29per, mname1 = "JHI-Hv50k-2016-272251")
plotPXG(HvDRR29per, "JHI-Hv50k-2016-272251")

find.marker(HvDRR29per, 2, 224.7)
effectplot(HvDRR29per, mname1 = "JHI-Hv50k-2016-126298")
summary(HvDRR29per, mname1 = "JHI-Hv50k-2016-126298")
plotPXG(HvDRR29per, "JHI-Hv50k-2016-126298")

summary(fitqtl(HvDRR29per, qtl = HvDRR29rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3 + Q4))

HvDRR29CI2 <- lodint(HvDRR29rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR29CI22 <- lodint(HvDRR29rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)
HvDRR29CI4 <- lodint(HvDRR29rqtl, chr = 4, qtl.index = 3, expandtomarkers = TRUE)
HvDRR29CI222 <- lodint(HvDRR29rqtl, chr = 2, qtl.index = 4, expandtomarkers = TRUE)

pt(-10.681, 106, lower.tail = FALSE)
#p<0.01
pt(-5.143, 106, lower.tail = FALSE)
#p<0.01
pt(5.135, 106, lower.tail = FALSE)
#p<0.01
pt(3.773, 106, lower.tail = FALSE)
#p<0.01


#HvDRR30####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
#library(qtl)
HvDRR30 <- Populations$HvDRR30

HvDRR30$geno$'1'$data<- apply(HvDRR30$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'2'$data<- apply(HvDRR30$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'3'$data<- apply(HvDRR30$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'4'$data<- apply(HvDRR30$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'5'$data<- apply(HvDRR30$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'6'$data<- apply(HvDRR30$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'7'$data<- apply(HvDRR30$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'1'$data<- apply(HvDRR30$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'2'$data<- apply(HvDRR30$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'3'$data<- apply(HvDRR30$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'4'$data<- apply(HvDRR30$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'5'$data<- apply(HvDRR30$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'6'$data<- apply(HvDRR30$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'7'$data<- apply(HvDRR30$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR30)[1] <- "riself"


summary(HvDRR30)
plot(HvDRR30)
HvDRR30 <- subset(HvDRR30,ind=c(1:nrow(HvDRR30$pheno))[!is.na(HvDRR30$pheno)])

HvDRR30per <- calc.genoprob(HvDRR30, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR30out <- scanone(HvDRR30per, method = "hk")
plot(HvDRR30out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR30permout <- scanone(HvDRR30per, method = "hk", n.perm=4000)
summary(HvDRR30out, perms = HvDRR30permout, alpha = 0.05, pvalues = TRUE)
HvDRR30qtl1 <- makeqtl(HvDRR30per, chr = c(2,7), pos = c(44,77), what = "prob")
HvDRR30outc27 <- addqtl(HvDRR30per, qtl = HvDRR30qtl1, method = "hk")
summary(HvDRR30outc27, perms = HvDRR30permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR30out, HvDRR30outc27, chr = -2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR30out, HvDRR30outc27, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)
HvDRR30qtl2 <- makeqtl(HvDRR30per, c(2,7,3,4), c(44,77,238,230), what = "prob")
HvDRR30outc2734 <- addqtl(HvDRR30per, qtl = HvDRR30qtl2, method = "hk")
summary(HvDRR30outc2734, perms = HvDRR30permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR30out, HvDRR30outc2734, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)
HvDRR30qtl3 <- makeqtl(HvDRR30per, c(2,7,3,4,2), c(44,77,238,230,123), what = "prob")
HvDRR30outc27342 <- addqtl(HvDRR30per, qtl = HvDRR30qtl3, method = "hk")
summary(HvDRR30outc27342, perms = HvDRR30permout, alpha = 0.05, pvalues = TRUE)

HvDRR30rqtl <- refineqtl(HvDRR30per, qtl = HvDRR30qtl3, method = "hk", verbose = FALSE)
summary(HvDRR30rqtl)
summary(fitqtl(HvDRR30per, qtl = HvDRR30rqtl, method = "hk"), pvalues = FALSE)


find.marker(HvDRR30per, 2, 43)
effectplot(HvDRR30per, mname1 = "JHI-Hv50k-2016-73566")
summary(HvDRR30per, mname1 = "JHI-Hv50k-2016-73566")
plotPXG(HvDRR30per, "JHI-Hv50k-2016-73566")

find.marker(HvDRR30per, 7, 76.8)
effectplot(HvDRR30per, mname1 = "SCRI_RS_160723")
summary(HvDRR30per, mname1 = "SCRI_RS_160723")
plotPXG(HvDRR30per, "SCRI_RS_160723")

find.marker(HvDRR30per, 3, 230.3)
effectplot(HvDRR30per, mname1 = "JHI-Hv50k-2016-212590")
summary(HvDRR30per, mname1 = "JHI-Hv50k-2016-212590")
plotPXG(HvDRR30per, "JHI-Hv50k-2016-212590")

find.marker(HvDRR30per, 4, 233)
effectplot(HvDRR30per, mname1 = "JHI-Hv50k-2016-272270")
summary(HvDRR30per, mname1 = "JHI-Hv50k-2016-272270")
plotPXG(HvDRR30per, "JHI-Hv50k-2016-272270")

find.marker(HvDRR30per, 2, 122.8)
effectplot(HvDRR30per, mname1 = "JHI-Hv50k-2016-97657")
summary(HvDRR30per, mname1 = "JHI-Hv50k-2016-97657")
plotPXG(HvDRR30per, "JHI-Hv50k-2016-97657")

summary(fitqtl(HvDRR30per, qtl = HvDRR30rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3 + Q4 + Q5))

HvDRR30CI2 <- lodint(HvDRR30rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR30CI7 <- lodint(HvDRR30rqtl, chr = 7, qtl.index = 2, expandtomarkers = TRUE)
HvDRR30CI3 <- lodint(HvDRR30rqtl, chr = 3, qtl.index = 3, expandtomarkers = TRUE)
HvDRR30CI4 <- lodint(HvDRR30rqtl, chr = 4, qtl.index = 4, expandtomarkers = TRUE)
HvDRR30CI22 <- lodint(HvDRR30rqtl, chr = 2, qtl.index = 5, expandtomarkers = TRUE)

pt(9.697, 117, lower.tail = FALSE)
#p<0.01
pt(6.521, 117, lower.tail = FALSE)
#p<0.01
pt(3.809, 117, lower.tail = FALSE)
#p<0.01
pt(-6.826, 117, lower.tail = FALSE)
#p<0.01
pt(4.615, 117, lower.tail = FALSE)
#p<0.01








#HvDRR31####
Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR31 <- Populations$HvDRR31

HvDRR31$geno$'1'$data<- apply(HvDRR31$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'2'$data<- apply(HvDRR31$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'3'$data<- apply(HvDRR31$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'4'$data<- apply(HvDRR31$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'5'$data<- apply(HvDRR31$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'6'$data<- apply(HvDRR31$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'7'$data<- apply(HvDRR31$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'1'$data<- apply(HvDRR31$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'2'$data<- apply(HvDRR31$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'3'$data<- apply(HvDRR31$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'4'$data<- apply(HvDRR31$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'5'$data<- apply(HvDRR31$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'6'$data<- apply(HvDRR31$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'7'$data<- apply(HvDRR31$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR31)[1] <- "riself"


summary(HvDRR31)
plot(HvDRR31)
HvDRR31 <- subset(HvDRR31,ind=c(1:nrow(HvDRR31$pheno))[!is.na(HvDRR31$pheno)])


HvDRR31per <- calc.genoprob(HvDRR31, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR31out <- scanone(HvDRR31per, method = "hk")
plot(HvDRR31out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR31permout <- scanone(HvDRR31per, method = "hk", n.perm=15000)
summary(HvDRR31out, perms = HvDRR31permout, alpha = 0.05, pvalues = TRUE)
HvDRR31qtl1 <- makeqtl(HvDRR31per, chr = 5, pos = 213, what = "prob")
HvDRR31outc5 <- addqtl(HvDRR31per, qtl = HvDRR31qtl1, method = "hk")
summary(HvDRR31outc5, perms = HvDRR31permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR31out, HvDRR31outc5, chr = -5, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR31out, HvDRR31outc5, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)
HvDRR31qtl2 <- makeqtl(HvDRR31per, chr = c(5,1,2,5), pos = c(213, 266, 175, 212), what = "prob")
HvDRR31outc5125 <- addqtl(HvDRR31per, qtl = HvDRR31qtl2, method = "hk")
summary(HvDRR31outc5125, perms = HvDRR31permout, alpha = 0.05, pvalues = TRUE)
HvDRR31qtl3 <- makeqtl(HvDRR31per, chr = c(5,1,2,5,5), pos = c(213, 266, 175, 212, 172), what = "prob")
HvDRR31outc51255 <- addqtl(HvDRR31per, qtl = HvDRR31qtl3, method = "hk")
summary(HvDRR31outc51255, perms = HvDRR31permout, alpha = 0.05, pvalues = TRUE)

HvDRR31rqtl <- refineqtl(HvDRR31per, qtl = HvDRR31qtl3, method = "hk", verbose = FALSE)
summary(HvDRR31rqtl)
summary(fitqtl(HvDRR31per, qtl = HvDRR31rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR31per, 5, 213)
effectplot(HvDRR31per, mname1 = "JHI-Hv50k-2016-338907")
summary(HvDRR31per, mname1 = "JHI-Hv50k-2016-338907")
plotPXG(HvDRR31per, "JHI-Hv50k-2016-338907")

find.marker(HvDRR31per, 1, 263.4)
effectplot(HvDRR31per, mname1 = "JHI-Hv50k-2016-55370")
summary(HvDRR31per, mname1 = "JHI-Hv50k-2016-55370")
plotPXG(HvDRR31per, "JHI-Hv50k-2016-55370")

find.marker(HvDRR31per, 2, 155)
effectplot(HvDRR31per, mname1 = "JHI-Hv50k-2016-91422")
summary(HvDRR31per, mname1 = "JHI-Hv50k-2016-91422")
plotPXG(HvDRR31per, "JHI-Hv50k-2016-91422")

find.marker(HvDRR31per, 5, 210)
effectplot(HvDRR31per, mname1 = "JHI-Hv50k-2016-338278")
summary(HvDRR31per, mname1 = "JHI-Hv50k-2016-338278")
plotPXG(HvDRR31per, "JHI-Hv50k-2016-338278")

find.marker(HvDRR31per, 5, 172)
effectplot(HvDRR31per, mname1 = "JHI-Hv50k-2016-321854")
summary(HvDRR31per, mname1 = "JHI-Hv50k-2016-321854")
plotPXG(HvDRR31per, "JHI-Hv50k-2016-321854")

summary(fitqtl(HvDRR31per, qtl = HvDRR31rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1 + Q2 + Q3 + Q4 + Q5))

HvDRR31CI5 <- lodint(HvDRR31rqtl, chr = 5, qtl.index = 1, expandtomarkers = TRUE)
HvDRR31CI1 <- lodint(HvDRR31rqtl, chr = 1, qtl.index = 2, expandtomarkers = TRUE)
HvDRR31CI2 <- lodint(HvDRR31rqtl, chr = 2, qtl.index = 3, expandtomarkers = TRUE)
HvDRR31CI55 <- lodint(HvDRR31rqtl, chr = 5, qtl.index = 4, expandtomarkers = TRUE)
HvDRR31CI555 <- lodint(HvDRR31rqtl, chr = 5, qtl.index = 5, expandtomarkers = TRUE)

pt(6.054, 123, lower.tail = FALSE)
#p<0.01
pt(5.502, 123, lower.tail = FALSE)
#p<0.01
pt(-5.075, 123, lower.tail = FALSE)
#p<0.01
pt(-4.797, 123, lower.tail = FALSE)
#p<0.01
pt(4.922, 123, lower.tail = FALSE)
#p<0.01





#HvDRR32####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR32 <- Populations$HvDRR32

HvDRR32$geno$'1'$data<- apply(HvDRR32$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'2'$data<- apply(HvDRR32$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'3'$data<- apply(HvDRR32$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'4'$data<- apply(HvDRR32$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'5'$data<- apply(HvDRR32$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'6'$data<- apply(HvDRR32$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'7'$data<- apply(HvDRR32$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'1'$data<- apply(HvDRR32$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'2'$data<- apply(HvDRR32$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'3'$data<- apply(HvDRR32$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'4'$data<- apply(HvDRR32$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'5'$data<- apply(HvDRR32$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'6'$data<- apply(HvDRR32$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'7'$data<- apply(HvDRR32$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR32)[1] <- "riself"

summary(HvDRR32)
plot(HvDRR32)
HvDRR32 <- subset(HvDRR32,ind=c(1:nrow(HvDRR32$pheno))[!is.na(HvDRR32$pheno)])

testmarkername <- as.data.frame(HvDRR32$geno$`4`$map)
which(row.names(testmarkername) == 'JHI-Hv50k-2016-271694')
which(row.names(testmarkername) == 'JHI-Hv50k-2016-271206')
which(row.names(testmarkername) == 'JHI-Hv50k-2016-272270')
marker2drop <- markernames(HvDRR32, chr = 4)[c(100,106,107)]
HvDRR32 <- drop.markers(HvDRR32, marker2drop)

HvDRR32per <- calc.genoprob(HvDRR32, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR32out <- scanone(HvDRR32per, method = "hk")
plot(HvDRR32out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR32permout <- scanone(HvDRR32per, method = "hk", n.perm=4000)
summary(HvDRR32out, perms = HvDRR32permout, alpha = 0.05, pvalues = TRUE)
HvDRR32qtl1 <- makeqtl(HvDRR32per, chr = 4, pos = 203, what = "prob")
HvDRR32outc4 <- addqtl(HvDRR32per, qtl = HvDRR32qtl1, method = "hk")
summary(HvDRR32outc4, perms = HvDRR32permout, alpha = 0.05, pvalues = TRUE)
HvDRR32qtl2 <- makeqtl(HvDRR32per, chr =  c(4,2), pos = c(203,101), what = "prob")
HvDRR32outc42 <- addqtl(HvDRR32per, qtl = HvDRR32qtl2, method = "hk")
summary(HvDRR32outc42, perms = HvDRR32permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR32out, HvDRR32outc42, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)

HvDRR32rqtl <- refineqtl(HvDRR32per, qtl = HvDRR32qtl2, method = "hk", verbose = FALSE)
summary(HvDRR32rqtl)
summary(fitqtl(HvDRR32per, qtl = HvDRR32rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR32per, 4, 202.5)
effectplot(HvDRR32per, mname1 = "JHI-Hv50k-2016-275319")
summary(HvDRR32per, mname1 = "JHI-Hv50k-2016-275319")
plotPXG(HvDRR32per, "JHI-Hv50k-2016-275319")

find.marker(HvDRR32per, 2, 101.3)
effectplot(HvDRR32per, mname1 = "JHI-Hv50k-2016-99627")
summary(HvDRR32per, mname1 = "JHI-Hv50k-2016-99627")
plotPXG(HvDRR32per, "JHI-Hv50k-2016-99627")

summary(fitqtl(HvDRR32per, qtl = HvDRR32rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2))

HvDRR32CI4 <- lodint(HvDRR32rqtl, chr = 4, qtl.index = 1, expandtomarkers = TRUE)
HvDRR32CI2 <- lodint(HvDRR32rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)

pt(-5.774, 35, lower.tail = FALSE)
#p<0.01
pt(4.339, 35, lower.tail = FALSE)
#p<0.01




#HvDRR33####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR33 <- Populations$HvDRR33

HvDRR33$geno$'1'$data<- apply(HvDRR33$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'2'$data<- apply(HvDRR33$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'3'$data<- apply(HvDRR33$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'4'$data<- apply(HvDRR33$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'5'$data<- apply(HvDRR33$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'6'$data<- apply(HvDRR33$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'7'$data<- apply(HvDRR33$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

HvDRR33$geno$'1'$data<- apply(HvDRR33$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'2'$data<- apply(HvDRR33$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'3'$data<- apply(HvDRR33$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'4'$data<- apply(HvDRR33$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'5'$data<- apply(HvDRR33$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'6'$data<- apply(HvDRR33$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'7'$data<- apply(HvDRR33$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR33)[1] <- "riself"

summary(HvDRR33)
plot(HvDRR33)
HvDRR33 <- subset(HvDRR33,ind=c(1:nrow(HvDRR33$pheno))[!is.na(HvDRR33$pheno)])

HvDRR33per <- calc.genoprob(HvDRR33, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR33out <- scanone(HvDRR33per, method = "hk")
plot(HvDRR33out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR33permout <- scanone(HvDRR33per, method = "hk", n.perm = 4000)
summary(HvDRR33out, perms = HvDRR33permout, alpha = 0.05, pvalues = TRUE)
HvDRR33qtl1 <- makeqtl(HvDRR33per, chr = 2, pos = 61.8, what = "prob")
HvDRR33outc2 <- addqtl(HvDRR33per, qtl = HvDRR33qtl1, method = "hk")
summary(HvDRR33outc2, perms = HvDRR33permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR33out, HvDRR33outc2, chr = -2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR33out, HvDRR33outc2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)

HvDRR33rqtl <- refineqtl(HvDRR33per, qtl = HvDRR33qtl1, method = "hk", verbose = FALSE)
summary(HvDRR33rqtl)
summary(fitqtl(HvDRR33per, qtl = HvDRR33rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR33per, 2, 61.8)
effectplot(HvDRR33per, mname1 = "JHI-Hv50k-2016-73510")
summary(HvDRR33per, mname1 = "JHI-Hv50k-2016-73510")
plotPXG(HvDRR33per, "JHI-Hv50k-2016-73510")

summary(fitqtl(HvDRR33per, qtl = HvDRR33rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR33CI2 <- lodint(HvDRR33rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(-5.66, 78, lower.tail = FALSE)
#p<0.01





#HvDRR34####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR34 <- Populations$HvDRR34

HvDRR34$geno$'1'$data<- apply(HvDRR34$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'2'$data<- apply(HvDRR34$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'3'$data<- apply(HvDRR34$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'4'$data<- apply(HvDRR34$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'5'$data<- apply(HvDRR34$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'6'$data<- apply(HvDRR34$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'7'$data<- apply(HvDRR34$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
#to eliminate heterozygous
HvDRR34$geno$'1'$data<- apply(HvDRR34$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'2'$data<- apply(HvDRR34$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'3'$data<- apply(HvDRR34$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'4'$data<- apply(HvDRR34$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'5'$data<- apply(HvDRR34$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'6'$data<- apply(HvDRR34$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'7'$data<- apply(HvDRR34$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR34)[1] <- "riself"


summary(HvDRR34)
plot(HvDRR34)
HvDRR34 <- subset(HvDRR34,ind=c(1:nrow(HvDRR34$pheno))[!is.na(HvDRR34$pheno)])

HvDRR34per <- calc.genoprob(HvDRR34, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR34out <- scanone(HvDRR34per, method = "hk")
plot(HvDRR34out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR34permout <- scanone(HvDRR34per, method = "hk", n.perm=4000)
summary(HvDRR34out, perms = HvDRR34permout, alpha = 0.05, pvalues = TRUE)
HvDRR34qtl1 <- makeqtl(HvDRR34per, chr = 2, pos = 289, what = "prob")
HvDRR34outc2 <- addqtl(HvDRR34per, qtl = HvDRR34qtl1, method = "hk")
summary(HvDRR34outc2, perms = HvDRR34permout, alpha = 0.05, pvalues = TRUE)

HvDRR34rqtl <- refineqtl(HvDRR34per, qtl = HvDRR34qtl1, method = "hk", verbose = FALSE)
summary(HvDRR34rqtl)
summary(fitqtl(HvDRR34per, qtl = HvDRR34rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR34per, 2, 289)
effectplot(HvDRR34per, mname1 = "JHI-Hv50k-2016-132317")
summary(HvDRR34per, mname1 = "JHI-Hv50k-2016-132317")
plotPXG(HvDRR34per, "JHI-Hv50k-2016-132317")

summary(fitqtl(HvDRR34per, qtl = HvDRR34rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR34CI2 <- lodint(HvDRR34rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(-4.739, 32, lower.tail = FALSE)
#p<0.01



#HVDRR35####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR35 <- Populations$HvDRR35

HvDRR35$geno$'1'$data<- apply(HvDRR35$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'2'$data<- apply(HvDRR35$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'3'$data<- apply(HvDRR35$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'4'$data<- apply(HvDRR35$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'5'$data<- apply(HvDRR35$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'6'$data<- apply(HvDRR35$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'7'$data<- apply(HvDRR35$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'1'$data<- apply(HvDRR35$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'2'$data<- apply(HvDRR35$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'3'$data<- apply(HvDRR35$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'4'$data<- apply(HvDRR35$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'5'$data<- apply(HvDRR35$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'6'$data<- apply(HvDRR35$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'7'$data<- apply(HvDRR35$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR35)[1] <- "riself"

summary(HvDRR35)
plot(HvDRR35)
HvDRR35 <- subset(HvDRR35,ind=c(1:nrow(HvDRR35$pheno))[!is.na(HvDRR35$pheno)])

HvDRR35per <- calc.genoprob(HvDRR35, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR35out <- scanone(HvDRR35per, method = "hk")
plot(HvDRR35out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR35permout <- scanone(HvDRR35per, method = "hk", n.perm=4000)
summary(HvDRR35out, perms = HvDRR35permout, alpha = 0.05, pvalues = TRUE)
HvDRR35qtl1 <- makeqtl(HvDRR35per, chr = 5, pos = 257, what = "prob")
HvDRR35outc5 <- addqtl(HvDRR35per, qtl = HvDRR35qtl1, method = "hk")
summary(HvDRR35outc5, perms = HvDRR35permout, alpha = 0.05, pvalues = TRUE)
HvDRR35qtl2 <- makeqtl(HvDRR35per, chr = c(5,2), pos = c(257,176), what = "prob")
HvDRR35outc52 <- addqtl(HvDRR35per, qtl = HvDRR35qtl2, method = "hk")
summary(HvDRR35outc52, perms = HvDRR35permout, alpha = 0.05, pvalues = TRUE)

HvDRR35rqtl <- refineqtl(HvDRR35per, qtl = HvDRR35qtl2, method = "hk", verbose = FALSE)
summary(HvDRR35rqtl)
summary(fitqtl(HvDRR35per, qtl = HvDRR35rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR35per, 5, 253.4)
effectplot(HvDRR35per, mname1 = "JHI-Hv50k-2016-333666")
summary(HvDRR35per, mname1 = "JHI-Hv50k-2016-333666")
plotPXG(HvDRR35per, "JHI-Hv50k-2016-333666")

find.marker(HvDRR35per, 2,175.8)
effectplot(HvDRR35per, mname1 = "JHI-Hv50k-2016-99962")
summary(HvDRR35per, mname1 = "JHI-Hv50k-2016-99962")
plotPXG(HvDRR35per, "JHI-Hv50k-2016-99962")

summary(fitqtl(HvDRR35per, qtl = HvDRR35rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1 + Q2))

HvDRR35CI5 <- lodint(HvDRR35rqtl, chr = 5, qtl.index = 1, expandtomarkers = TRUE)
HvDRR35CI2 <- lodint(HvDRR35rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)

pt(5.109, 82, lower.tail = FALSE)
#p<0.01
pt(-4.181, 82, lower.tail = FALSE)
#not significant

#HvDRR36####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR36 <- Populations$HvDRR36

HvDRR36$geno$'1'$data<- apply(HvDRR36$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'2'$data<- apply(HvDRR36$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'3'$data<- apply(HvDRR36$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'4'$data<- apply(HvDRR36$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'5'$data<- apply(HvDRR36$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'6'$data<- apply(HvDRR36$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'7'$data<- apply(HvDRR36$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
#to eliminate heterozygous
HvDRR36$geno$'1'$data<- apply(HvDRR36$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'2'$data<- apply(HvDRR36$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'3'$data<- apply(HvDRR36$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'4'$data<- apply(HvDRR36$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'5'$data<- apply(HvDRR36$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'6'$data<- apply(HvDRR36$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'7'$data<- apply(HvDRR36$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR36)[1] <- "riself"

summary(HvDRR36)
plot(HvDRR36)
HvDRR36 <- subset(HvDRR36,ind=c(1:nrow(HvDRR36$pheno))[!is.na(HvDRR36$pheno)])

testmarkername <- as.data.frame(HvDRR36$geno$`4`$map)
which(row.names(testmarkername) == 'JHI-Hv50k-2016-276443')
which(row.names(testmarkername) == 'JHI-Hv50k-2016-274011')
which(row.names(testmarkername) == 'SCRI_RS_119390')
marker2drop <- markernames(HvDRR36, chr = 4)[c(137,138,139)]
HvDRR36 <- drop.markers(HvDRR36, marker2drop)

HvDRR36per <- calc.genoprob(HvDRR36, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR36out <- scanone(HvDRR36per, method = "hk")
plot(HvDRR36out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR36permout <- scanone(HvDRR36per, method = "hk", n.perm=4000)
summary(HvDRR36out, perms = HvDRR36permout, alpha = 0.05, pvalues = TRUE)
HvDRR36qtl1 <- makeqtl(HvDRR36per, chr = 4, pos = 197, what = "prob")
HvDRR36outc4 <- addqtl(HvDRR36per, qtl = HvDRR36qtl1, method = "hk")
summary(HvDRR36outc4, perms = HvDRR36permout, alpha = 0.05, pvalues = TRUE)

HvDRR36rqtl <- refineqtl(HvDRR36per, qtl = HvDRR36qtl1, method = "hk", verbose = FALSE)
summary(HvDRR36rqtl)
summary(fitqtl(HvDRR36per, qtl = HvDRR36rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR36per, qtl = HvDRR36rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))


find.marker(HvDRR36per, 4, 197)
effectplot(HvDRR36per, mname1 = "JHI-Hv50k-2016-273434")
summary(HvDRR36per, mname1 = "JHI-Hv50k-2016-273434")
plotPXG(HvDRR36per, "JHI-Hv50k-2016-273434")

HvDRR36CI4 <- lodint(HvDRR36rqtl, chr = 4, qtl.index = 1, expandtomarkers = TRUE)

pt(-5.221, 56, lower.tail = FALSE)
#p<0.01


#HvDRR37####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR37 <- Populations$HvDRR37

HvDRR37$geno$'1'$data<- apply(HvDRR37$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'2'$data<- apply(HvDRR37$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'3'$data<- apply(HvDRR37$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'4'$data<- apply(HvDRR37$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'5'$data<- apply(HvDRR37$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'6'$data<- apply(HvDRR37$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'7'$data<- apply(HvDRR37$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'1'$data<- apply(HvDRR37$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'2'$data<- apply(HvDRR37$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'3'$data<- apply(HvDRR37$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'4'$data<- apply(HvDRR37$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'5'$data<- apply(HvDRR37$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'6'$data<- apply(HvDRR37$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'7'$data<- apply(HvDRR37$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR37)[1] <- "riself"

summary(HvDRR37)
plot(HvDRR37)
HvDRR37 <- subset(HvDRR37,ind=c(1:nrow(HvDRR37$pheno))[!is.na(HvDRR37$pheno)])

HvDRR37per <- calc.genoprob(HvDRR37, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR37out <- scanone(HvDRR37per, method = "hk")
plot(HvDRR37out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR37permout <- scanone(HvDRR37per, method = "hk", n.perm=4000)
summary(HvDRR37out, perms = HvDRR37permout, alpha = 0.05, pvalues = TRUE)
HvDRR37qtl1 <- makeqtl(HvDRR37per, chr = 2, pos = 75.9, what = "prob")
HvDRR37outc2 <- addqtl(HvDRR37per, qtl = HvDRR37qtl1, method = "hk")
summary(HvDRR37outc2, perms = HvDRR37permout, alpha = 0.05, pvalues = TRUE)

HvDRR37rqtl <- refineqtl(HvDRR37per, qtl = HvDRR37qtl1, method = "hk", verbose = FALSE)
summary(HvDRR37rqtl)
summary(fitqtl(HvDRR37per, qtl = HvDRR37rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR37per, qtl = HvDRR37rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))


find.marker(HvDRR37per, 2, 75.9)
effectplot(HvDRR37per, mname1 = "SCRI_RS_207244")
summary(HvDRR37per, mname1 = "SCRI_RS_207244")
plotPXG(HvDRR37per, "SCRI_RS_207244")

HvDRR37CI2 <- lodint(HvDRR37rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(-4.074, 60, lower.tail = FALSE)
#p<0.01


#HvDRR38####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR38 <- Populations$HvDRR38

HvDRR38$geno$'1'$data<- apply(HvDRR38$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'2'$data<- apply(HvDRR38$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'3'$data<- apply(HvDRR38$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'4'$data<- apply(HvDRR38$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'5'$data<- apply(HvDRR38$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'6'$data<- apply(HvDRR38$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'7'$data<- apply(HvDRR38$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
#to eliminate heterozygous
HvDRR38$geno$'1'$data<- apply(HvDRR38$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'2'$data<- apply(HvDRR38$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'3'$data<- apply(HvDRR38$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'4'$data<- apply(HvDRR38$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'5'$data<- apply(HvDRR38$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'6'$data<- apply(HvDRR38$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'7'$data<- apply(HvDRR38$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR38)[1] <- "riself"

summary(HvDRR38)
plot(HvDRR38)
HvDRR38 <- subset(HvDRR38,ind=c(1:nrow(HvDRR38$pheno))[!is.na(HvDRR38$pheno)])


HvDRR38per <- calc.genoprob(HvDRR38, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR38out <- scanone(HvDRR38per, method = "hk")
plot(HvDRR38out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR38permout <- scanone(HvDRR38per, method = "hk", n.perm=4000)
summary(HvDRR38out, perms = HvDRR38permout, alpha = 0.05, pvalues = TRUE)
HvDRR38qtl1 <- makeqtl(HvDRR38per, chr = 2, pos = 360, what = "prob")
HvDRR38outc2 <- addqtl(HvDRR38per, qtl = HvDRR38qtl1, method = "hk")
summary(HvDRR38outc2, perms = HvDRR38permout, alpha = 0.05, pvalues = TRUE)

HvDRR38rqtl <- refineqtl(HvDRR38per, qtl = HvDRR38qtl1, method = "hk", verbose = FALSE)
summary(HvDRR38rqtl)
summary(fitqtl(HvDRR38per, qtl = HvDRR38rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR38per, qtl = HvDRR38rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))


find.marker(HvDRR38per, 2, 360)
effectplot(HvDRR38per, mname1 = "SCRI_RS_147230")
summary(HvDRR38per, mname1 = "SCRI_RS_147230")
plotPXG(HvDRR38per, "SCRI_RS_147230")

HvDRR38CI2 <- lodint(HvDRR38rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(-4.645, 74, lower.tail = FALSE)
#p<0.01

#HvDRR39####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
#library(qtl)
HvDRR39 <- Populations$HvDRR39

HvDRR39$geno$'1'$data<- apply(HvDRR39$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'2'$data<- apply(HvDRR39$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'3'$data<- apply(HvDRR39$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'4'$data<- apply(HvDRR39$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'5'$data<- apply(HvDRR39$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'6'$data<- apply(HvDRR39$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'7'$data<- apply(HvDRR39$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'1'$data<- apply(HvDRR39$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'2'$data<- apply(HvDRR39$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'3'$data<- apply(HvDRR39$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'4'$data<- apply(HvDRR39$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'5'$data<- apply(HvDRR39$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'6'$data<- apply(HvDRR39$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'7'$data<- apply(HvDRR39$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR39)[1] <- "riself"

summary(HvDRR39)
plot(HvDRR39)
HvDRR39 <- subset(HvDRR39,ind=c(1:nrow(HvDRR39$pheno))[!is.na(HvDRR39$pheno)])

HvDRR39per <- calc.genoprob(HvDRR39, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR39out <- scanone(HvDRR39per, method = "hk")
plot(HvDRR39out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR39permout <- scanone(HvDRR39per, method = "hk", n.perm=15000)
summary(HvDRR39out, perms = HvDRR39permout, alpha = 0.05, pvalues = TRUE)
HvDRR39qtl1 <- makeqtl(HvDRR39per, chr = 2, pos = 240, what = "prob")
HvDRR39outc2 <- addqtl(HvDRR39per, qtl = HvDRR39qtl1, method = "hk")
summary(HvDRR39outc2, perms = HvDRR39permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR39out, HvDRR39outc2, chr = -2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR39out, HvDRR39outc2, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)
HvDRR39qtl2 <- makeqtl(HvDRR39per, c(2,2), c(240, 41.6), what = "prob")
HvDRR39outc22 <- addqtl(HvDRR39per, qtl = HvDRR39qtl2, method = "hk")
summary(HvDRR39outc22, perms = HvDRR39permout, alpha = 0.05, pvalues = TRUE)
HvDRR39qtl3 <- makeqtl(HvDRR39per, c(2,2,2), c(240, 41.6, 121), what = "prob")
HvDRR39outc222 <- addqtl(HvDRR39per, qtl = HvDRR39qtl3, method = "hk")
summary(HvDRR39outc222, perms = HvDRR39permout, alpha = 0.05, pvalues = TRUE)

HvDRR39rqtl <- refineqtl(HvDRR39per, qtl = HvDRR39qtl3, method = "hk", verbose = FALSE)
summary(HvDRR39rqtl)
summary(fitqtl(HvDRR39per, qtl = HvDRR39rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR39per, 2, 239.7)
effectplot(HvDRR39per, mname1 = "SCRI_RS_121952")
summary(HvDRR39per, mname1 = "SCRI_RS_121952")
plotPXG(HvDRR39per, "SCRI_RS_121952")

find.marker(HvDRR39per, 2, 23)
effectplot(HvDRR39per, mname1 = "JHI-Hv50k-2016-69230")
summary(HvDRR39per, mname1 = "JHI-Hv50k-2016-69230")
plotPXG(HvDRR39per, "JHI-Hv50k-2016-69230")

find.marker(HvDRR39per, 2, 121)
effectplot(HvDRR39per, mname1 = "JHI-Hv50k-2016-94897")
summary(HvDRR39per, mname1 = "JHI-Hv50k-2016-94897")
plotPXG(HvDRR39per, "JHI-Hv50k-2016-94897")

summary(fitqtl(HvDRR39per, qtl = HvDRR39rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3))

HvDRR39CI2 <- lodint(HvDRR39rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR39CI22 <- lodint(HvDRR39rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)
HvDRR39CI222 <- lodint(HvDRR39rqtl, chr = 2, qtl.index = 3, expandtomarkers = TRUE)

pt(-6.185, 101, lower.tail = FALSE)
#p<0.01
pt(5.819, 101, lower.tail = FALSE)
#p<0.01
pt(-4.925, 101, lower.tail = FALSE)
#p<0.01


#HvDRR40####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR40 <- Populations$HvDRR40

HvDRR40$geno$'1'$data<- apply(HvDRR40$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'2'$data<- apply(HvDRR40$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'3'$data<- apply(HvDRR40$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'4'$data<- apply(HvDRR40$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'5'$data<- apply(HvDRR40$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'6'$data<- apply(HvDRR40$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'7'$data<- apply(HvDRR40$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'1'$data<- apply(HvDRR40$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'2'$data<- apply(HvDRR40$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'3'$data<- apply(HvDRR40$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'4'$data<- apply(HvDRR40$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'5'$data<- apply(HvDRR40$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'6'$data<- apply(HvDRR40$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'7'$data<- apply(HvDRR40$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
class(HvDRR40)[1] <- "riself"

summary(HvDRR40)
plot(HvDRR40)
HvDRR40 <- subset(HvDRR40,ind=c(1:nrow(HvDRR40$pheno))[!is.na(HvDRR40$pheno)])

HvDRR40per <- calc.genoprob(HvDRR40, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR40out <- scanone(HvDRR40per, method = "hk")
plot(HvDRR40out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR40permout <- scanone(HvDRR40per, method = "hk", n.perm=4000)
summary(HvDRR40out, perms = HvDRR40permout, alpha = 0.05, pvalues = TRUE)


#HvDRR41####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR41 <- Populations$HvDRR41

HvDRR41$geno$'1'$data<- apply(HvDRR41$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'2'$data<- apply(HvDRR41$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'3'$data<- apply(HvDRR41$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'4'$data<- apply(HvDRR41$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'5'$data<- apply(HvDRR41$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'6'$data<- apply(HvDRR41$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'7'$data<- apply(HvDRR41$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'1'$data<- apply(HvDRR41$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'2'$data<- apply(HvDRR41$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'3'$data<- apply(HvDRR41$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'4'$data<- apply(HvDRR41$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'5'$data<- apply(HvDRR41$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'6'$data<- apply(HvDRR41$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'7'$data<- apply(HvDRR41$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR41)[1] <- "riself"

summary(HvDRR41)
plot(HvDRR41)
HvDRR41 <- subset(HvDRR41,ind=c(1:nrow(HvDRR41$pheno))[!is.na(HvDRR41$pheno)])

HvDRR41per <- calc.genoprob(HvDRR41, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR41out <- scanone(HvDRR41per, method = "hk")
plot(HvDRR41out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR41permout <- scanone(HvDRR41per, method = "hk", n.perm=4000)
summary(HvDRR41out, perms = HvDRR41permout, alpha = 0.05, pvalues = TRUE)
HvDRR41qtl1 <- makeqtl(HvDRR41per, chr = 2, pos = 45.8, what = "prob")
HvDRR41outc2 <- addqtl(HvDRR41per, qtl = HvDRR41qtl1, method = "hk")
summary(HvDRR41outc2, perms = HvDRR41permout, alpha = 0.05, pvalues = TRUE)
HvDRR41qtl2 <- makeqtl(HvDRR41per, chr = c(2,2), pos = c(45.8,130), what = "prob")
HvDRR41outc22 <- addqtl(HvDRR41per, qtl = HvDRR41qtl2, method = "hk")
summary(HvDRR41outc22, perms = HvDRR41permout, alpha = 0.05, pvalues = TRUE)
HvDRR41qtl3 <- makeqtl(HvDRR41per, chr = c(2,2,7), pos = c(45.8,130,52), what = "prob")
HvDRR41outc227 <- addqtl(HvDRR41per, qtl = HvDRR41qtl3, method = "hk")
summary(HvDRR41outc227, perms = HvDRR41permout, alpha = 0.05, pvalues = TRUE)

HvDRR41rqtl <- refineqtl(HvDRR41per, qtl = HvDRR41qtl3, method = "hk", verbose = FALSE)
summary(HvDRR41rqtl)
summary(fitqtl(HvDRR41per, qtl = HvDRR41rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR41per, 2, 45.8)
effectplot(HvDRR41per, mname1 = "JHI-Hv50k-2016-73417")
summary(HvDRR41per, mname1 = "JHI-Hv50k-2016-73417")
plotPXG(HvDRR41per, "JHI-Hv50k-2016-73417")

find.marker(HvDRR41per, 2, 129.6)
effectplot(HvDRR41per, mname1 = "BOPA2_12_11278")
summary(HvDRR41per, mname1 = "BOPA2_12_11278")
plotPXG(HvDRR41per, "BOPA2_12_11278")

find.marker(HvDRR41per, 7, 52)
effectplot(HvDRR41per, mname1 = "JHI-Hv50k-2016-460797")
summary(HvDRR41per, mname1 = "JHI-Hv50k-2016-460797")
plotPXG(HvDRR41per, "JHI-Hv50k-2016-460797")

summary(fitqtl(HvDRR41per, qtl = HvDRR41rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1 + Q2 + Q3))

HvDRR41CI2 <- lodint(HvDRR41rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR41CI22 <- lodint(HvDRR41rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)
HvDRR41CI7 <- lodint(HvDRR41rqtl, chr = 7, qtl.index = 3, expandtomarkers = TRUE)
pt(8.660, 89, lower.tail = FALSE)
pt(5.999, 89, lower.tail = FALSE)
pt(5.259, 89, lower.tail = FALSE)



#HvDRR42####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR42 <- Populations$HvDRR42

HvDRR42$geno$'1'$data<- apply(HvDRR42$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'2'$data<- apply(HvDRR42$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'3'$data<- apply(HvDRR42$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'4'$data<- apply(HvDRR42$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'5'$data<- apply(HvDRR42$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'6'$data<- apply(HvDRR42$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'7'$data<- apply(HvDRR42$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'1'$data<- apply(HvDRR42$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'2'$data<- apply(HvDRR42$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'3'$data<- apply(HvDRR42$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'4'$data<- apply(HvDRR42$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'5'$data<- apply(HvDRR42$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'6'$data<- apply(HvDRR42$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'7'$data<- apply(HvDRR42$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR42)[1] <- "riself"

summary(HvDRR42)
plot(HvDRR42)
HvDRR42 <- subset(HvDRR42,ind=c(1:nrow(HvDRR42$pheno))[!is.na(HvDRR42$pheno)])

HvDRR42per <- calc.genoprob(HvDRR42, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR42out <- scanone(HvDRR42per, method = "hk")
plot(HvDRR42out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR42permout <- scanone(HvDRR42per, method = "hk", n.perm=4000)
summary(HvDRR42out, perms = HvDRR42permout, alpha = 0.05, pvalues = TRUE)
HvDRR42qtl1 <- makeqtl(HvDRR42per, chr = c(2,3,4), pos = c(68.3,291,155), what = "prob")
HvDRR42outc234 <- addqtl(HvDRR42per, qtl = HvDRR42qtl1, method = "hk")
summary(HvDRR42outc234, perms = HvDRR42permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR42out, HvDRR42outc234, chr = c(-3,-4), col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR42out, HvDRR42outc234, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)
HvDRR42qtl2 <- makeqtl(HvDRR42per, chr = c(2,3,4,5), pos = c(68.3,291,155,456), what = "prob")
HvDRR42outc2345 <- addqtl(HvDRR42per, qtl = HvDRR42qtl2, method = "hk")
summary(HvDRR42outc2345, perms = HvDRR42permout, alpha = 0.05, pvalues = TRUE)

HvDRR42rqtl <- refineqtl(HvDRR42per, qtl = HvDRR42qtl2, method = "hk", verbose = FALSE)
summary(HvDRR42rqtl)
summary(fitqtl(HvDRR42per, qtl = HvDRR42rqtl, method = "hk"), pvalues = FALSE)
plot(HvDRR42rqtl)

find.marker(HvDRR42per, 2, 67.6)
effectplot(HvDRR42per, mname1 = "JHI-Hv50k-2016-72896")
summary(HvDRR42per, mname1 = "JHI-Hv50k-2016-72896")
plotPXG(HvDRR42per, "JHI-Hv50k-2016-72896")

find.marker(HvDRR42per, 3, 286.4)
effectplot(HvDRR42per, mname1 = "JHI-Hv50k-2016-196425")
summary(HvDRR42per, mname1 = "JHI-Hv50k-2016-196425")
plotPXG(HvDRR42per, "JHI-Hv50k-2016-196425")

find.marker(HvDRR42per, 4, 155.1)
effectplot(HvDRR42per, mname1 = "BOPA2_12_30455")
summary(HvDRR42per, mname1 = "BOPA2_12_30455")
plotPXG(HvDRR42per, "BOPA2_12_30455")

find.marker(HvDRR42per, 5, 456)
effectplot(HvDRR42per, mname1 = "JHI-Hv50k-2016-338482")
summary(HvDRR42per, mname1 = "JHI-Hv50k-2016-338482")
plotPXG(HvDRR42per, "JHI-Hv50k-2016-338482")

summary(fitqtl(HvDRR42per, qtl = HvDRR42rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3 + Q4))

HvDRR42CI2 <- lodint(HvDRR42rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR42CI3 <- lodint(HvDRR42rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)
HvDRR42CI4 <- lodint(HvDRR42rqtl, chr = 4, qtl.index = 3, expandtomarkers = TRUE)
HvDRR42CI5 <- lodint(HvDRR42rqtl, chr = 5, qtl.index = 4, expandtomarkers = TRUE)

pt(4.430, 67, lower.tail = FALSE)
#p<0.01
pt(-6.340, 67, lower.tail = FALSE)
#not significant
pt(-3.563, 67, lower.tail = FALSE)
#p<0.01
pt(-5.516, 67, lower.tail = FALSE)
#p<0.01




#HvDRR43####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR43 <- Populations$HvDRR43

HvDRR43$geno$'1'$data<- apply(HvDRR43$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'2'$data<- apply(HvDRR43$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'3'$data<- apply(HvDRR43$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'4'$data<- apply(HvDRR43$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'5'$data<- apply(HvDRR43$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'6'$data<- apply(HvDRR43$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'7'$data<- apply(HvDRR43$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'1'$data<- apply(HvDRR43$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'2'$data<- apply(HvDRR43$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'3'$data<- apply(HvDRR43$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'4'$data<- apply(HvDRR43$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'5'$data<- apply(HvDRR43$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'6'$data<- apply(HvDRR43$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'7'$data<- apply(HvDRR43$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR43)[1] <- "riself"

summary(HvDRR43)
plot(HvDRR43)
HvDRR43 <- subset(HvDRR43,ind=c(1:nrow(HvDRR43$pheno))[!is.na(HvDRR43$pheno)])

HvDRR43per <- calc.genoprob(HvDRR43, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR43out <- scanone(HvDRR43per, method = "hk")
plot(HvDRR43out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR43permout <- scanone(HvDRR43per, method = "hk", n.perm=4000)
summary(HvDRR43out, perms = HvDRR43permout, alpha = 0.05, pvalues = TRUE)
HvDRR43qtl1 <- makeqtl(HvDRR43per, chr = 4, pos = 131, what = "prob")
HvDRR43outc4 <- addqtl(HvDRR43per, qtl = HvDRR43qtl1, method = "hk")
summary(HvDRR43outc4, perms = HvDRR43permout, alpha = 0.05, pvalues = TRUE)

HvDRR43rqtl <- refineqtl(HvDRR43per, qtl = HvDRR43qtl1, method = "hk", verbose = FALSE)
summary(HvDRR43rqtl)
summary(fitqtl(HvDRR43per, qtl = HvDRR43rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR43per, 4, 131.1)
effectplot(HvDRR43per, mname1 = "JHI-Hv50k-2016-268870")
summary(HvDRR43per, mname1 = "JHI-Hv50k-2016-268870")
plotPXG(HvDRR43per, "JHI-Hv50k-2016-268870")

summary(fitqtl(HvDRR43per, qtl = HvDRR43rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR43CI <- lodint(HvDRR43rqtl, chr = 4, qtl.index = 1, expandtomarkers = TRUE)

pt(-8.281, 119, lower.tail = FALSE)
#p<0.01



#HvDRR44####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR44 <- Populations$HvDRR44

HvDRR44$geno$'1'$data<- apply(HvDRR44$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'2'$data<- apply(HvDRR44$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'3'$data<- apply(HvDRR44$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'4'$data<- apply(HvDRR44$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'5'$data<- apply(HvDRR44$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'6'$data<- apply(HvDRR44$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'7'$data<- apply(HvDRR44$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'1'$data<- apply(HvDRR44$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'2'$data<- apply(HvDRR44$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'3'$data<- apply(HvDRR44$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'4'$data<- apply(HvDRR44$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'5'$data<- apply(HvDRR44$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'6'$data<- apply(HvDRR44$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'7'$data<- apply(HvDRR44$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR44)[1] <- "riself"

summary(HvDRR44)
plot(HvDRR44)
HvDRR44 <- subset(HvDRR44,ind=c(1:nrow(HvDRR44$pheno))[!is.na(HvDRR44$pheno)])

HvDRR44per <- calc.genoprob(HvDRR44, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR44out <- scanone(HvDRR44per, method = "hk")
plot(HvDRR44out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR44permout <- scanone(HvDRR44per, method = "hk", n.perm = 4000)
summary(HvDRR44out, perms = HvDRR44permout, alpha = 0.05, pvalues = TRUE)
HvDRR44qtl1 <- makeqtl(HvDRR44per, chr = 5, pos = 250, what = "prob")
HvDRR44outc5 <- addqtl(HvDRR44per, qtl = HvDRR44qtl1, method = "hk")
summary(HvDRR44outc5, perms = HvDRR44permout, alpha = 0.05, pvalues = TRUE)
HvDRR44qtl2 <- makeqtl(HvDRR44per, chr = c(5,7), pos = c(250,74.7), what = "prob")
HvDRR44outc57 <- addqtl(HvDRR44per, qtl = HvDRR44qtl2, method = "hk")
summary(HvDRR44outc57, perms = HvDRR44permout, alpha = 0.05, pvalues = TRUE)

HvDRR44rqtl <- refineqtl(HvDRR44per, qtl = HvDRR44qtl2, method = "hk", verbose = FALSE)
summary(HvDRR44rqtl)
summary(fitqtl(HvDRR44per, qtl = HvDRR44rqtl, method = "hk"), pvalues = FALSE)
summary(fitqtl(HvDRR44per, qtl = HvDRR44rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1+Q2))

find.marker(HvDRR44per, 5, 261)
effectplot(HvDRR44per, mname1 = "JHI-Hv50k-2016-334016")
summary(HvDRR44per, mname1 = "JHI-Hv50k-2016-334016")
plotPXG(HvDRR44per, "JHI-Hv50k-2016-334016")

find.marker(HvDRR44per, 7, 74.7)
effectplot(HvDRR44per, mname1 = "JHI-Hv50k-2016-459818")
summary(HvDRR44per, mname1 = "JHI-Hv50k-2016-459818")
plotPXG(HvDRR44per, "JHI-Hv50k-2016-459818")

HvDRR44CI5 <- lodint(HvDRR44rqtl, chr = 5, qtl.index = 1, expandtomarkers = TRUE)
HvDRR44CI7 <- lodint(HvDRR44rqtl, chr = 7, qtl.index = 2, expandtomarkers = TRUE)

pt(5.383, 87, lower.tail = FALSE)
pt(4.728, 87, lower.tail = FALSE)


#HvDRR45####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR45 <- Populations$HvDRR45

HvDRR45$geno$'1'$data<- apply(HvDRR45$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'2'$data<- apply(HvDRR45$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'3'$data<- apply(HvDRR45$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'4'$data<- apply(HvDRR45$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'5'$data<- apply(HvDRR45$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'6'$data<- apply(HvDRR45$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'7'$data<- apply(HvDRR45$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'1'$data<- apply(HvDRR45$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'2'$data<- apply(HvDRR45$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'3'$data<- apply(HvDRR45$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'4'$data<- apply(HvDRR45$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'5'$data<- apply(HvDRR45$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'6'$data<- apply(HvDRR45$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'7'$data<- apply(HvDRR45$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR45)[1] <- "riself"

summary(HvDRR45)
plot(HvDRR45)
HvDRR45 <- subset(HvDRR45,ind=c(1:nrow(HvDRR45$pheno))[!is.na(HvDRR45$pheno)])

HvDRR45per <- calc.genoprob(HvDRR45, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR45out <- scanone(HvDRR45per, method = "hk")
plot(HvDRR45out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR45permout <- scanone(HvDRR45per, method = "hk", n.perm=4000)
summary(HvDRR45out, perms = HvDRR45permout, alpha = 0.05, pvalues = TRUE)
HvDRR45qtl1 <- makeqtl(HvDRR45per, chr = 5, pos = 200, what = "prob")
HvDRR45outc2 <- addqtl(HvDRR45per, qtl = HvDRR45qtl1, method = "hk")
summary(HvDRR45outc2, perms = HvDRR45permout, alpha = 0.05, pvalues = TRUE)

HvDRR45rqtl <- refineqtl(HvDRR45per, qtl = HvDRR45qtl1, method = "hk", verbose = FALSE)
summary(HvDRR45rqtl)
summary(fitqtl(HvDRR45per, qtl = HvDRR45rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR45per, 5, 200)
effectplot(HvDRR45per, mname1 = "JHI-Hv50k-2016-337093")
summary(HvDRR45per, mname1 = "JHI-Hv50k-2016-337093")
plotPXG(HvDRR45per, "JHI-Hv50k-2016-337093")

summary(fitqtl(HvDRR45per, qtl = HvDRR45rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1))

HvDRR45CI5 <- lodint(HvDRR45rqtl, chr = 5, qtl.index = 1, expandtomarkers = TRUE)

pt(5.731, 75, lower.tail = FALSE)
#p<0.01

#HvDRR46####


Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR46 <- Populations$HvDRR46

HvDRR46$geno$'1'$data<- apply(HvDRR46$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'2'$data<- apply(HvDRR46$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'3'$data<- apply(HvDRR46$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'4'$data<- apply(HvDRR46$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'5'$data<- apply(HvDRR46$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'6'$data<- apply(HvDRR46$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'7'$data<- apply(HvDRR46$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'1'$data<- apply(HvDRR46$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'2'$data<- apply(HvDRR46$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'3'$data<- apply(HvDRR46$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'4'$data<- apply(HvDRR46$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'5'$data<- apply(HvDRR46$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'6'$data<- apply(HvDRR46$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'7'$data<- apply(HvDRR46$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR46)[1] <- "riself"


summary(HvDRR46)
plot(HvDRR46)
HvDRR46 <- subset(HvDRR46,ind=c(1:nrow(HvDRR46$pheno))[!is.na(HvDRR46$pheno)])

HvDRR46per <- calc.genoprob(HvDRR46, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR46out <- scanone(HvDRR46per, method = "hk")
plot(HvDRR46out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR46permout <- scanone(HvDRR46per, method = "hk", n.perm=4000)
summary(HvDRR46out, perms = HvDRR46permout, alpha = 0.05, pvalues = TRUE)
HvDRR46qtl1 <- makeqtl(HvDRR46per, chr = 5, pos = 355, what = "prob")
HvDRR46outc5 <- addqtl(HvDRR46per, qtl = HvDRR46qtl1, method = "hk")
summary(HvDRR46outc5, perms = HvDRR46permout, alpha = 0.05, pvalues = TRUE)

HvDRR46rqtl <- refineqtl(HvDRR46per, qtl = HvDRR46qtl1, method = "hk", verbose = FALSE)
summary(HvDRR46rqtl)
summary(fitqtl(HvDRR46per, qtl = HvDRR46rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR46per, 5, 355.3)
effectplot(HvDRR46per, mname1 = "JHI-Hv50k-2016-335345")
summary(HvDRR46per, mname1 = "JHI-Hv50k-2016-335345")
plotPXG(HvDRR46per, "JHI-Hv50k-2016-335345")

summary(fitqtl(HvDRR46per, qtl = HvDRR46rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1))

HvDRR46CI5 <- lodint(HvDRR46rqtl, chr = 5, qtl.index = 1, expandtomarkers = TRUE)

pt(4.097, 32, lower.tail = FALSE)
#p<0.01

#HvDRR47####


Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR47 <- Populations$HvDRR47

HvDRR47$geno$'1'$data<- apply(HvDRR47$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'2'$data<- apply(HvDRR47$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'3'$data<- apply(HvDRR47$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'4'$data<- apply(HvDRR47$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'5'$data<- apply(HvDRR47$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'6'$data<- apply(HvDRR47$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'7'$data<- apply(HvDRR47$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'1'$data<- apply(HvDRR47$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'2'$data<- apply(HvDRR47$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'3'$data<- apply(HvDRR47$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'4'$data<- apply(HvDRR47$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'5'$data<- apply(HvDRR47$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'6'$data<- apply(HvDRR47$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'7'$data<- apply(HvDRR47$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR47)[1] <- "riself"


summary(HvDRR47)
plot(HvDRR47)
HvDRR47 <- subset(HvDRR47,ind=c(1:nrow(HvDRR47$pheno))[!is.na(HvDRR47$pheno)])

testmarkername <- as.data.frame(HvDRR47$geno$`4`$map)
which(row.names(testmarkername) == 'JHI-Hv50k-2016-275849')
which(row.names(testmarkername) == 'JHI-Hv50k-2016-276097')
marker2drop <- markernames(HvDRR47, chr = 4)[c(394:395,402:412)]
HvDRR47 <- drop.markers(HvDRR47, marker2drop)

HvDRR47per <- calc.genoprob(HvDRR47, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR47out <- scanone(HvDRR47per, method = "hk")
plot(HvDRR47out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR47permout <- scanone(HvDRR47per, method = "hk", n.perm = 10000)
summary(HvDRR47out, perms = HvDRR47permout, alpha = 0.05, pvalues = TRUE)

HvDRR47qtl1 <- makeqtl(HvDRR47per, c(2,4,5,7), c(40,240,266,344), what = "prob")
HvDRR47outc245 <- addqtl(HvDRR47per, qtl = HvDRR47qtl1, method = "hk")
summary(HvDRR47outc245, perms = HvDRR47permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR47out, HvDRR47outc245, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)

HvDRR47rqtl <- refineqtl(HvDRR47per, qtl = HvDRR47qtl1, method = "hk", verbose = FALSE)
summary(HvDRR47rqtl)
summary(fitqtl(HvDRR47per, qtl = HvDRR47rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR47per, 2, 40.7)
effectplot(HvDRR47per, mname1 = "JHI-Hv50k-2016-73422")
summary(HvDRR47per, mname1 = "JHI-Hv50k-2016-73422")
plotPXG(HvDRR47per, "JHI-Hv50k-2016-73422")

find.marker(HvDRR47per, 4, 242)
effectplot(HvDRR47per, mname1 = "SCRI_RS_221876")
summary(HvDRR47per, mname1 = "SCRI_RS_221876")
plotPXG(HvDRR47per, "SCRI_RS_221876")

find.marker(HvDRR47per, 5, 273)
effectplot(HvDRR47per, mname1 = "SCRI_RS_193529")
summary(HvDRR47per, mname1 = "SCRI_RS_193529")
plotPXG(HvDRR47per, "SCRI_RS_193529")

find.marker(HvDRR47per, 7, 200)
effectplot(HvDRR47per, mname1 = "JHI-Hv50k-2016-475698")
summary(HvDRR47per, mname1 = "JHI-Hv50k-2016-475698")
plotPXG(HvDRR47per, "JHI-Hv50k-2016-475698")

summary(fitqtl(HvDRR47per, qtl = HvDRR47rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1 + Q2 + Q3 + Q4))

HvDRR47CI2 <- lodint(HvDRR47rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR47CI4 <- lodint(HvDRR47rqtl, chr = 4, qtl.index = 2, expandtomarkers = TRUE)
HvDRR47CI5 <- lodint(HvDRR47rqtl, chr = 5, qtl.index = 3, expandtomarkers = TRUE)
HvDRR47CI7 <- lodint(HvDRR47rqtl, chr = 7, qtl.index = 4, expandtomarkers = TRUE)

pt(8.266, 95, lower.tail = FALSE)
#p<0.01
pt(-7.702, 95, lower.tail = FALSE)
#p<0.01
pt(-7-100, 95, lower.tail = FALSE)
#p<0.01
pt(1605, 95, lower.tail = FALSE)
#p<0.01

#HvDRR48####

Populations <- readRDS("Population_analysis_Flowering_Time.RDS")
library(qtl)
HvDRR48 <- Populations$HvDRR48

HvDRR48$geno$'1'$data<- apply(HvDRR48$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'2'$data<- apply(HvDRR48$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'3'$data<- apply(HvDRR48$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'4'$data<- apply(HvDRR48$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'5'$data<- apply(HvDRR48$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'6'$data<- apply(HvDRR48$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'7'$data<- apply(HvDRR48$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'1'$data<- apply(HvDRR48$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'2'$data<- apply(HvDRR48$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'3'$data<- apply(HvDRR48$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'4'$data<- apply(HvDRR48$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'5'$data<- apply(HvDRR48$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'6'$data<- apply(HvDRR48$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'7'$data<- apply(HvDRR48$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR48)[1] <- "riself"

summary(HvDRR48)
plot(HvDRR48)
HvDRR48 <- subset(HvDRR48,ind=c(1:nrow(HvDRR48$pheno))[!is.na(HvDRR48$pheno)])

testmarkername <- as.data.frame(HvDRR48$geno$`2`$map)
which(row.names(testmarkername) == 'JHI-Hv50k-2016-72531')
marker2drop <- markernames(HvDRR48, chr = 2)[46]
HvDRR48 <- drop.markers(HvDRR48, marker2drop)
HvDRR48per <- calc.genoprob(HvDRR48, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR48out <- scanone(HvDRR48per, method = "hk")
plot(HvDRR48out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR48permout <- scanone(HvDRR48per, method = "hk", n.perm=4000)
summary(HvDRR48out, perms = HvDRR48permout, alpha = 0.05, pvalues = TRUE)
HvDRR48qtl1 <- makeqtl(HvDRR48per, chr = 2, pos = 50.6, what = "prob")
HvDRR48outc2 <- addqtl(HvDRR48per, qtl = HvDRR48qtl1, method = "hk")
summary(HvDRR48outc2, perms = HvDRR48permout, alpha = 0.05, pvalues = TRUE)

HvDRR48rqtl <- refineqtl(HvDRR48per, qtl = HvDRR48qtl1, method = "hk", verbose = FALSE)
summary(HvDRR48rqtl)
summary(fitqtl(HvDRR48per, qtl = HvDRR48rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR48per, 2, 50.6)
effectplot(HvDRR48per, mname1 = "JHI-Hv50k-2016-75253")
summary(HvDRR48per, mname1 = "JHI-Hv50k-2016-75253")
plotPXG(HvDRR48per, "JHI-Hv50k-2016-75253")

summary(fitqtl(HvDRR48per, qtl = HvDRR48rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1))

HvDRR48CI2 <- lodint(HvDRR48rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(-4.635, 82, lower.tail = FALSE)
#p<0.01


#PLANT HEIGHT##################
setwd("/QTLanalysis/Plant_height/")

#HvDRR02####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR02 <- Populations$HvDRR02

HvDRR02$geno$'1'$data<- apply(HvDRR02$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'2'$data<- apply(HvDRR02$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'3'$data<- apply(HvDRR02$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'4'$data<- apply(HvDRR02$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'5'$data<- apply(HvDRR02$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'6'$data<- apply(HvDRR02$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'7'$data<- apply(HvDRR02$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR02$geno$'1'$data<- apply(HvDRR02$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'2'$data<- apply(HvDRR02$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'3'$data<- apply(HvDRR02$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'4'$data<- apply(HvDRR02$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'5'$data<- apply(HvDRR02$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'6'$data<- apply(HvDRR02$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR02$geno$'7'$data<- apply(HvDRR02$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR02)[1] <- "riself"

summary(HvDRR02)
plot(HvDRR02)
HvDRR02 <- subset(HvDRR02,ind=c(1:nrow(HvDRR02$pheno))[!is.na(HvDRR02$pheno)])

HvDRR02per <- calc.genoprob(HvDRR02, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR02out <- scanone(HvDRR02per, method = "hk")
plot(HvDRR02out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR02permout <- scanone(HvDRR02per, method = "hk", n.perm=4000)
summary(HvDRR02out, perms = HvDRR02permout, alpha = 0.05, pvalues = TRUE)
HvDRR02qtl1 <- makeqtl(HvDRR02per, chr = c(2,5), pos = c(237, 214), what = "prob")
HvDRR02outc25 <- addqtl(HvDRR02per, qtl = HvDRR02qtl1, method = "hk")
summary(HvDRR02outc25, perms = HvDRR02permout, alpha = 0.05, pvalues = TRUE)
HvDRR02qtl2 <- makeqtl(HvDRR02per, chr = c(2,5,3), pos = c(237, 214, 196), what = "prob")
HvDRR02outc253 <- addqtl(HvDRR02per, qtl = HvDRR02qtl2, method = "hk")
summary(HvDRR02outc253, perms = HvDRR02permout, alpha = 0.05, pvalues = TRUE)
HvDRR02qtl3 <- makeqtl(HvDRR02per, chr = c(2,5,3,4), pos = c(237, 214, 196, 44.1), what = "prob")
HvDRR02outc2534 <- addqtl(HvDRR02per, qtl = HvDRR02qtl3, method = "hk")
summary(HvDRR02outc2534, perms = HvDRR02permout, alpha = 0.05, pvalues = TRUE)

HvDRR02rqtl <- refineqtl(HvDRR02per, qtl = HvDRR02qtl3, method = "hk", verbose = FALSE)
summary(HvDRR02rqtl)
summary(fitqtl(HvDRR02per, qtl = HvDRR02rqtl, method = "hk"), pvalues = FALSE)


find.marker(HvDRR02per, 2, 242.7)
effectplot(HvDRR02per, mname1 = "JHI-Hv50k-2016-132327")
summary(HvDRR02per, mname1 = "JHI-Hv50k-2016-132327")
plotPXG(HvDRR02per, "JHI-Hv50k-2016-132327")

find.marker(HvDRR02per, 5, 216.9)
effectplot(HvDRR02per, mname1 = "JHI-Hv50k-2016-338436")
summary(HvDRR02per, mname1 = "JHI-Hv50k-2016-338436")
plotPXG(HvDRR02per, "JHI-Hv50k-2016-338436")

find.marker(HvDRR02per, 3, 196.5)
effectplot(HvDRR02per, mname1 = "JHI-Hv50k-2016-207933")
summary(HvDRR02per, mname1 = "JHI-Hv50k-2016-207933")
plotPXG(HvDRR02per, "JHI-Hv50k-2016-207933")

find.marker(HvDRR02per, 4, 44.1)
effectplot(HvDRR02per, mname1 = "JHI-Hv50k-2016-231100")
summary(HvDRR02per, mname1 = "JHI-Hv50k-2016-231100")
plotPXG(HvDRR02per, "JHI-Hv50k-2016-231100")

summary(fitqtl(HvDRR02per, qtl = HvDRR02rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3 + Q4))

HvDRR02CI2 <- lodint(HvDRR02rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR02CI5 <- lodint(HvDRR02rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)
HvDRR02CI3 <- lodint(HvDRR02rqtl, chr = 3, qtl.index = 3, expandtomarkers = TRUE)
HvDRR02CI4 <- lodint(HvDRR02rqtl, chr = 4, qtl.index = 4, expandtomarkers = TRUE)

pt(-4.317, 133, lower.tail = FALSE)
#p<0.01
pt(-6.120, 133, lower.tail = FALSE)
#p<0.01
pt(4.745, 133, lower.tail = FALSE)
#p<0.01
pt(4.088, 133, lower.tail = FALSE)
#p<0.01

#HvDRR03####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR03 <- Populations$HvDRR03

HvDRR03$geno$'1'$data<- apply(HvDRR03$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'2'$data<- apply(HvDRR03$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'3'$data<- apply(HvDRR03$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'4'$data<- apply(HvDRR03$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'5'$data<- apply(HvDRR03$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'6'$data<- apply(HvDRR03$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'7'$data<- apply(HvDRR03$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR03$geno$'1'$data<- apply(HvDRR03$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'2'$data<- apply(HvDRR03$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'3'$data<- apply(HvDRR03$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'4'$data<- apply(HvDRR03$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'5'$data<- apply(HvDRR03$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'6'$data<- apply(HvDRR03$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR03$geno$'7'$data<- apply(HvDRR03$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR03)[1] <- "riself"

summary(HvDRR03)
plot(HvDRR03)
HvDRR03 <- subset(HvDRR03,ind=c(1:nrow(HvDRR03$pheno))[!is.na(HvDRR03$pheno)])

HvDRR03per <- calc.genoprob(HvDRR03, step = 1,  error.prob = 0.01, map.function = "haldane")
HvDRR03out <- scanone(HvDRR03per, method = "hk")
plot(HvDRR03out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR03permout <- scanone(HvDRR03per, method = "hk", n.perm = 4000)
summary(HvDRR03out, perms = HvDRR03permout, alpha = 0.05, pvalues = TRUE)

#HvDRR04####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR04 <- Populations$HvDRR04

HvDRR04$geno$'1'$data<- apply(HvDRR04$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR04$geno$'2'$data<- apply(HvDRR04$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR04$geno$'3'$data<- apply(HvDRR04$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR04$geno$'4'$data<- apply(HvDRR04$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR04$geno$'5'$data<- apply(HvDRR04$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR04$geno$'6'$data<- apply(HvDRR04$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR04$geno$'7'$data<- apply(HvDRR04$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR04$geno$'1'$data<- apply(HvDRR04$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR04$geno$'2'$data<- apply(HvDRR04$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR04$geno$'3'$data<- apply(HvDRR04$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR04$geno$'4'$data<- apply(HvDRR04$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR04$geno$'5'$data<- apply(HvDRR04$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR04$geno$'6'$data<- apply(HvDRR04$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR04$geno$'7'$data<- apply(HvDRR04$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR04)[1] <- "riself"


summary(HvDRR04)
plot(HvDRR04)
HvDRR04 <- subset(HvDRR04,ind=c(1:nrow(HvDRR04$pheno))[!is.na(HvDRR04$pheno)])

HvDRR04per <- calc.genoprob(HvDRR04, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR04out <- scanone(HvDRR04per, method = "hk")
plot(HvDRR04out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR04permout <- scanone(HvDRR04per, method = "hk", n.perm = 4000)
summary(HvDRR04out, perms = HvDRR04permout, alpha = 0.05, pvalues = TRUE)
HvDRR04qtl1 <- makeqtl(HvDRR04per, chr = 2, pos = 210, what = "prob")
HvDRR04outc2 <- addqtl(HvDRR04per, qtl = HvDRR04qtl1, method = "hk")
summary(HvDRR04outc2, perms = HvDRR04permout, alpha = 0.05, pvalues = TRUE)

HvDRR04rqtl <- refineqtl(HvDRR04per, qtl = HvDRR04qtl1, method = "hk", verbose = FALSE)
summary(HvDRR04rqtl)
summary(fitqtl(HvDRR04per, qtl = HvDRR04rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR04per, 2, 209.6)
effectplot(HvDRR04per, mname1 = "JHI-Hv50k-2016-116511")
summary(HvDRR04per, mname1 = "JHI-Hv50k-2016-116511")
plotPXG(HvDRR04per, "JHI-Hv50k-2016-116511")

summary(fitqtl(HvDRR04per, qtl = HvDRR04rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR04CI2 <- lodint(HvDRR04rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(-4.002, 140, lower.tail = FALSE)
#p<0.01
#HvDRR07####


Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR07 <- Populations$HvDRR07

HvDRR07$geno$'1'$data<- apply(HvDRR07$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'2'$data<- apply(HvDRR07$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'3'$data<- apply(HvDRR07$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'4'$data<- apply(HvDRR07$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'5'$data<- apply(HvDRR07$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'6'$data<- apply(HvDRR07$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'7'$data<- apply(HvDRR07$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR07$geno$'1'$data<- apply(HvDRR07$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'2'$data<- apply(HvDRR07$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'3'$data<- apply(HvDRR07$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'4'$data<- apply(HvDRR07$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'5'$data<- apply(HvDRR07$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'6'$data<- apply(HvDRR07$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR07$geno$'7'$data<- apply(HvDRR07$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR07)[1] <- "riself"

summary(HvDRR07)
plot(HvDRR07)
HvDRR07 <- subset(HvDRR07,ind=c(1:nrow(HvDRR07$pheno))[!is.na(HvDRR07$pheno)])

HvDRR07per <- calc.genoprob(HvDRR07, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR07out <- scanone(HvDRR07per, method = "hk")
plot(HvDRR07out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR07permout <- scanone(HvDRR07per, method = "hk", n.perm = 4000)
summary(HvDRR07out, perms = HvDRR07permout, alpha = 0.05, pvalues = TRUE)
HvDRR07qtl1 <- makeqtl(HvDRR07per, chr = c(2,3,4,7), pos = c(45.1,157.8,101.3,127.7), what = "prob")
HvDRR07outc2347 <- addqtl(HvDRR07per, qtl = HvDRR07qtl1, method = "hk")
summary(HvDRR07outc2347, perms = HvDRR07permout, alpha = 0.05, pvalues = TRUE)

HvDRR07rqtl <- refineqtl(HvDRR07per, qtl = HvDRR07qtl1, method = "hk", verbose = FALSE)
summary(HvDRR07rqtl)
summary(fitqtl(HvDRR07per, qtl = HvDRR07rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR07per, 2, 46)
effectplot(HvDRR07per, mname1 = "SCRI_RS_170337")
summary(HvDRR07per, mname1 = "SCRI_RS_170337")
plotPXG(HvDRR07per, "SCRI_RS_170337")

find.marker(HvDRR07per, 3, 158)
effectplot(HvDRR07per, mname1 = "JHI-Hv50k-2016-205562")
summary(HvDRR07per, mname1 = "JHI-Hv50k-2016-205562")
plotPXG(HvDRR07per, "JHI-Hv50k-2016-205562")

find.marker(HvDRR07per, 4, 124)
effectplot(HvDRR07per, mname1 = "JHI-Hv50k-2016-268537")
summary(HvDRR07per, mname1 = "JHI-Hv50k-2016-268537")
plotPXG(HvDRR07per, "JHI-Hv50k-2016-268537")

find.marker(HvDRR07per, 7, 129)
effectplot(HvDRR07per, mname1 = "JHI-Hv50k-2016-475138")
summary(HvDRR07per, mname1 = "JHI-Hv50k-2016-475138")
plotPXG(HvDRR07per, "JHI-Hv50k-2016-475138")

summary(fitqtl(HvDRR07per, qtl = HvDRR07rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3 + Q4))

HvDRR07CI2 <- lodint(HvDRR07rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR07CI3 <- lodint(HvDRR07rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)
HvDRR07CI4 <- lodint(HvDRR07rqtl, chr = 4, qtl.index = 3, expandtomarkers = TRUE)
HvDRR07CI7 <- lodint(HvDRR07rqtl, chr = 7, qtl.index = 4, expandtomarkers = TRUE)

pt(-7.987, 93, lower.tail = FALSE)
#p<0.01
pt(8.026, 93, lower.tail = FALSE)
#p<0.01
pt(-3.772, 93, lower.tail = FALSE)
#p<0.01
pt(5.647, 93, lower.tail = FALSE)
#p<0.01

#HvDRR08####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR08 <- Populations$HvDRR08


HvDRR08$geno$'1'$data<- apply(HvDRR08$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'2'$data<- apply(HvDRR08$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'3'$data<- apply(HvDRR08$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'4'$data<- apply(HvDRR08$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'5'$data<- apply(HvDRR08$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'6'$data<- apply(HvDRR08$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'7'$data<- apply(HvDRR08$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR08$geno$'1'$data<- apply(HvDRR08$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'2'$data<- apply(HvDRR08$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'3'$data<- apply(HvDRR08$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'4'$data<- apply(HvDRR08$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'5'$data<- apply(HvDRR08$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'6'$data<- apply(HvDRR08$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR08$geno$'7'$data<- apply(HvDRR08$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR08)[1] <- "riself"

summary(HvDRR08)
plot(HvDRR08)
HvDRR08 <- subset(HvDRR08,ind=c(1:nrow(HvDRR08$pheno))[!is.na(HvDRR08$pheno)])

testmarkername <- as.data.frame(HvDRR08$geno$`3`$map)
which(row.names(testmarkername) == 'BOPA1_5260-462')
marker2drop <- markernames(HvDRR08, chr = 3)[314]
HvDRR08 <- drop.markers(HvDRR08, marker2drop)
HvDRR08per <- calc.genoprob(HvDRR08, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR08out <- scanone(HvDRR08per, method = "hk")
plot(HvDRR08out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR08permout <- scanone(HvDRR08per, method = "hk", n.perm = 4000)
summary(HvDRR08out, perms = HvDRR08permout, alpha = 0.05, pvalues = TRUE)
HvDRR08qtl1 <- makeqtl(HvDRR08per, chr = c(2,3,4), pos = c(52.5,224.5,138), what = "prob")
HvDRR08outc234 <- addqtl(HvDRR08per, qtl = HvDRR08qtl1, method = "hk")
summary(HvDRR08outc234, perms = HvDRR08permout, alpha = 0.05, pvalues = TRUE)
HvDRR08qtl2 <- makeqtl(HvDRR08per, chr = c(2,3,4,5), pos = c(52.5,224.5,138,118), what = "prob")
HvDRR08outc2345 <- addqtl(HvDRR08per, qtl = HvDRR08qtl2, method = "hk")
summary(HvDRR08outc2345, perms = HvDRR08permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR08out, HvDRR08outc2345, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)

HvDRR08rqtl <- refineqtl(HvDRR08per, qtl = HvDRR08qtl2, method = "hk", verbose = FALSE)
summary(HvDRR08rqtl)
summary(fitqtl(HvDRR08per, qtl = HvDRR08rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR08per, 2, 57.6)
effectplot(HvDRR08per, mname1 = "JHI-Hv50k-2016-72853")
summary(HvDRR08per, mname1 = "JHI-Hv50k-2016-72853")
plotPXG(HvDRR08per, "JHI-Hv50k-2016-72853")

find.marker(HvDRR08per, 3, 223.6)
effectplot(HvDRR08per, mname1 = "JHI-Hv50k-2016-205259")
summary(HvDRR08per, mname1 = "JHI-Hv50k-2016-205259")
plotPXG(HvDRR08per, "JHI-Hv50k-2016-205259")

find.marker(HvDRR08per, 4, 135)
effectplot(HvDRR08per, mname1 = "SCRI_RS_176669")
summary(HvDRR08per, mname1 = "SCRI_RS_176669")
plotPXG(HvDRR08per, "SCRI_RS_176669")

find.marker(HvDRR08per, 5, 118.1)
effectplot(HvDRR08per, mname1 = "JHI-Hv50k-2016-307515")
summary(HvDRR08per, mname1 = "JHI-Hv50k-2016-307515")
plotPXG(HvDRR08per, "JHI-Hv50k-2016-307515")

HvDRR08CI2 <- lodint(HvDRR08rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR08CI3 <- lodint(HvDRR08rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)
HvDRR08CI4 <- lodint(HvDRR08rqtl, chr = 4, qtl.index = 3, expandtomarkers = TRUE)
HvDRR08CI5 <- lodint(HvDRR08rqtl, chr = 5, qtl.index = 4, expandtomarkers = TRUE)

summary(fitqtl(HvDRR08per, qtl = HvDRR08rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3 + Q4))

pt(6.462, 90, lower.tail = FALSE)
#p<0.01
pt(-6.566, 90, lower.tail = FALSE)
#not significant
pt(-5.944, 90, lower.tail = FALSE)
#p<0.01
pt(4.453, 90, lower.tail = FALSE)
#p<0.01

#HvDRR09####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR09 <- Populations$HvDRR09

HvDRR09$geno$'1'$data<- apply(HvDRR09$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'2'$data<- apply(HvDRR09$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'3'$data<- apply(HvDRR09$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'4'$data<- apply(HvDRR09$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'5'$data<- apply(HvDRR09$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'6'$data<- apply(HvDRR09$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'7'$data<- apply(HvDRR09$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR09$geno$'1'$data<- apply(HvDRR09$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'2'$data<- apply(HvDRR09$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'3'$data<- apply(HvDRR09$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'4'$data<- apply(HvDRR09$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'5'$data<- apply(HvDRR09$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'6'$data<- apply(HvDRR09$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR09$geno$'7'$data<- apply(HvDRR09$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
class(HvDRR09)[1] <- "riself"

summary(HvDRR09)
plot(HvDRR09)
HvDRR09 <- subset(HvDRR09,ind=c(1:nrow(HvDRR09$pheno))[!is.na(HvDRR09$pheno)])

HvDRR09per <- calc.genoprob(HvDRR09, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR09out <- scanone(HvDRR09per, method = "hk")
plot(HvDRR09out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR09permout <- scanone(HvDRR09per, method = "hk", n.perm = 4000)
summary(HvDRR09out, perms = HvDRR09permout, alpha = 0.05, pvalues = TRUE)

#HvDRR10####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR10 <- Populations$HvDRR10

HvDRR10$geno$'1'$data<- apply(HvDRR10$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'2'$data<- apply(HvDRR10$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'3'$data<- apply(HvDRR10$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'4'$data<- apply(HvDRR10$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'5'$data<- apply(HvDRR10$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'6'$data<- apply(HvDRR10$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR10$geno$'7'$data<- apply(HvDRR10$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
#to eliminate heterozygous
HvDRR10$geno$'1'$data<- apply(HvDRR10$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'2'$data<- apply(HvDRR10$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'3'$data<- apply(HvDRR10$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'4'$data<- apply(HvDRR10$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'5'$data<- apply(HvDRR10$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'6'$data<- apply(HvDRR10$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR10$geno$'7'$data<- apply(HvDRR10$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR10)[1] <- "riself"


summary(HvDRR10)
plot(HvDRR10)
HvDRR10 <- subset(HvDRR10,ind=c(1:nrow(HvDRR10$pheno))[!is.na(HvDRR10$pheno)])

HvDRR10per <- calc.genoprob(HvDRR10, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR10out <- scanone(HvDRR10per, method = "hk")
plot(HvDRR10out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR10permout <- scanone(HvDRR10per, method = "hk", n.perm = 10000)
summary(HvDRR10out, perms = HvDRR10permout, alpha = 0.05, pvalues = TRUE)
HvDRR10qtl1 <- makeqtl(HvDRR10per, chr = 3, pos = 124, what = "prob")
HvDRR10outc3 <- addqtl(HvDRR10per, qtl = HvDRR10qtl1, method = "hk")
summary(HvDRR10outc3, perms = HvDRR10permout, alpha = 0.05, pvalues = TRUE)
HvDRR10qtl2 <- makeqtl(HvDRR10per, chr = c(3,2), pos = c(124,189), what = "prob")
HvDRR10outc32 <- addqtl(HvDRR10per, qtl = HvDRR10qtl2, method = "hk")
summary(HvDRR10outc32, perms = HvDRR10permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR10out, HvDRR10outc32, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)

HvDRR10rqtl <- refineqtl(HvDRR10per, qtl = HvDRR10qtl2, method = "hk", verbose = FALSE)
summary(HvDRR10rqtl)
summary(fitqtl(HvDRR10per, qtl = HvDRR10rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR10per, 3, 123.7)
effectplot(HvDRR10per, mname1 = "JHI-Hv50k-2016-205404")
summary(HvDRR10per, mname1 = "JHI-Hv50k-2016-205404")
plotPXG(HvDRR10per, "JHI-Hv50k-2016-205404")

find.marker(HvDRR10per, 2, 189.2)
effectplot(HvDRR10per, mname1 = "JHI-Hv50k-2016-130123")
summary(HvDRR10per, mname1 = "JHI-Hv50k-2016-130123")
plotPXG(HvDRR10per, "JHI-Hv50k-2016-130123")

summary(fitqtl(HvDRR10per, qtl = HvDRR10rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2))

HvDRR10CI3 <- lodint(HvDRR10rqtl, chr = 3, qtl.index = 1, expandtomarkers = TRUE)
HvDRR10CI2 <- lodint(HvDRR10rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)

pt(-14.172, 86, lower.tail = FALSE)
#p<0.01
pt(4.713, 86, lower.tail = FALSE)
#p<0.01

#HvDRR11####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR11 <- Populations$HvDRR11

HvDRR11$geno$'1'$data<- apply(HvDRR11$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'2'$data<- apply(HvDRR11$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'3'$data<- apply(HvDRR11$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'4'$data<- apply(HvDRR11$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'5'$data<- apply(HvDRR11$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'6'$data<- apply(HvDRR11$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'7'$data<- apply(HvDRR11$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR11$geno$'1'$data<- apply(HvDRR11$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'2'$data<- apply(HvDRR11$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'3'$data<- apply(HvDRR11$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'4'$data<- apply(HvDRR11$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'5'$data<- apply(HvDRR11$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'6'$data<- apply(HvDRR11$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR11$geno$'7'$data<- apply(HvDRR11$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR11)[1] <- "riself"

summary(HvDRR11)
plot(HvDRR11)
HvDRR11 <- subset(HvDRR11,ind=c(1:nrow(HvDRR11$pheno))[!is.na(HvDRR11$pheno)])

testmarkername <- as.data.frame(HvDRR11$geno$`7`$map)
which(row.names(testmarkername) == 'JHI-Hv50k-2016-511608')
marker2drop <- markernames(HvDRR11, chr = 7)[201]
HvDRR11 <- drop.markers(HvDRR11, marker2drop)

testmarkername <- as.data.frame(HvDRR11$geno$`3`$map)
which(row.names(testmarkername) == 'SCRI_RS_150063')
marker2drop <- markernames(HvDRR11, chr = 3)[125]
HvDRR11 <- drop.markers(HvDRR11, marker2drop)

HvDRR11per <- calc.genoprob(HvDRR11, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR11out <- scanone(HvDRR11per, method = "hk")
plot(HvDRR11out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR11permout <- scanone(HvDRR11per, method = "hk", n.perm = 4000)
summary(HvDRR11out, perms = HvDRR11permout, alpha = 0.05, pvalues = TRUE)
HvDRR11qtl1 <- makeqtl(HvDRR11per, chr = c(2, 3), pos = c(166,112), what = "prob")
HvDRR11outc23 <- addqtl(HvDRR11per, qtl = HvDRR11qtl1, method = "hk")
summary(HvDRR11outc23, perms = HvDRR11permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR11out, HvDRR11outc23, chr = c(-2, -3), col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR11out, HvDRR11outc23, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)
HvDRR11qtl2 <- makeqtl(HvDRR11per, chr = c(2, 3, 2, 7), pos = c(166,112,137,77), what = "prob")
HvDRR11outc2327 <- addqtl(HvDRR11per, qtl = HvDRR11qtl2, method = "hk")
summary(HvDRR11outc2327, perms = HvDRR11permout, alpha = 0.05, pvalues = TRUE)
HvDRR11qtl3 <- makeqtl(HvDRR11per, chr = c(2, 3, 2, 7, 7), pos = c(166,112,137,77, 148), what = "prob")
HvDRR11outc23277 <- addqtl(HvDRR11per, qtl = HvDRR11qtl3, method = "hk")
summary(HvDRR11outc23277, perms = HvDRR11permout, alpha = 0.05, pvalues = TRUE)

HvDRR11rqtl <- refineqtl(HvDRR11per, qtl = HvDRR11qtl3, method = "hk", verbose = FALSE)
summary(HvDRR11rqtl)
summary(fitqtl(HvDRR11per, qtl = HvDRR11rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR11per, 2, 175.3)
effectplot(HvDRR11per, mname1 = "JHI-Hv50k-2016-135139")
summary(HvDRR11per, mname1 = "JHI-Hv50k-2016-135139")
plotPXG(HvDRR11per, "JHI-Hv50k-2016-135139")

find.marker(HvDRR11per, 3, 111)
effectplot(HvDRR11per, mname1 = "JHI-Hv50k-2016-205774")
summary(HvDRR11per, mname1 = "JHI-Hv50k-2016-205774")
plotPXG(HvDRR11per, "JHI-Hv50k-2016-205774")

find.marker(HvDRR11per, 2, 136.9)
effectplot(HvDRR11per, mname1 = "BOPA1_4833-420")
summary(HvDRR11per, mname1 = "BOPA1_4833-420")
plotPXG(HvDRR11per, "BOPA1_4833-420")

find.marker(HvDRR11per, 7, 76.8)
effectplot(HvDRR11per, mname1 = "JHI-Hv50k-2016-464308")
summary(HvDRR11per, mname1 = "JHI-Hv50k-2016-464308")
plotPXG(HvDRR11per, "JHI-Hv50k-2016-464308")

find.marker(HvDRR11per, 7, 160)
effectplot(HvDRR11per, mname1 = "JHI-Hv50k-2016-505814")
summary(HvDRR11per, mname1 = "JHI-Hv50k-2016-505814")
plotPXG(HvDRR11per, "JHI-Hv50k-2016-505814")

HvDRR11CI2 <- lodint(HvDRR11rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR11CI3 <- lodint(HvDRR11rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)
HvDRR11CI22 <- lodint(HvDRR11rqtl, chr = 2, qtl.index = 3, expandtomarkers = TRUE)
HvDRR11CI7 <- lodint(HvDRR11rqtl, chr = 7, qtl.index = 4, expandtomarkers = TRUE)
HvDRR11CI77 <- lodint(HvDRR11rqtl, chr = 7, qtl.index = 5, expandtomarkers = TRUE)

summary(fitqtl(HvDRR11per, qtl = HvDRR11rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3 + Q4 + Q5))

pt(7.819, 91, lower.tail = FALSE)
#p<0.01
pt(-11.031, 91, lower.tail = FALSE)
#p<0.01
pt(-2.697, 91, lower.tail = FALSE)
#p<0.01
pt(-4.215, 91, lower.tail = FALSE)
#p<0.01
pt(-4.253, 91, lower.tail = FALSE)
#p<0.01
#HvDRR12####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR12 <- Populations$HvDRR12

HvDRR12$geno$'1'$data<- apply(HvDRR12$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'2'$data<- apply(HvDRR12$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'3'$data<- apply(HvDRR12$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'4'$data<- apply(HvDRR12$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'5'$data<- apply(HvDRR12$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'6'$data<- apply(HvDRR12$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'7'$data<- apply(HvDRR12$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR12$geno$'1'$data<- apply(HvDRR12$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'2'$data<- apply(HvDRR12$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'3'$data<- apply(HvDRR12$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'4'$data<- apply(HvDRR12$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'5'$data<- apply(HvDRR12$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'6'$data<- apply(HvDRR12$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR12$geno$'7'$data<- apply(HvDRR12$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR12)[1] <- "riself"

summary(HvDRR12)
plot(HvDRR12)
HvDRR12 <- subset(HvDRR12,ind=c(1:nrow(HvDRR12$pheno))[!is.na(HvDRR12$pheno)])
HvDRR12per <- calc.genoprob(HvDRR12, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR12out <- scanone(HvDRR12per, method = "hk")
plot(HvDRR12out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR12permout <- scanone(HvDRR12per, method = "hk", n.perm = 4000)
summary(HvDRR12out, perms = HvDRR12permout, alpha = 0.05, pvalues = TRUE)

#HvDRR13####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR13 <- Populations$HvDRR13

HvDRR13$geno$'1'$data<- apply(HvDRR13$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'2'$data<- apply(HvDRR13$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'3'$data<- apply(HvDRR13$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'4'$data<- apply(HvDRR13$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'5'$data<- apply(HvDRR13$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'6'$data<- apply(HvDRR13$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'7'$data<- apply(HvDRR13$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR13$geno$'1'$data<- apply(HvDRR13$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'2'$data<- apply(HvDRR13$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'3'$data<- apply(HvDRR13$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'4'$data<- apply(HvDRR13$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'5'$data<- apply(HvDRR13$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'6'$data<- apply(HvDRR13$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR13$geno$'7'$data<- apply(HvDRR13$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR13)[1] <- "riself"


summary(HvDRR13)
plot(HvDRR13)
HvDRR13 <- subset(HvDRR13,ind=c(1:nrow(HvDRR13$pheno))[!is.na(HvDRR13$pheno)])
HvDRR13per <- calc.genoprob(HvDRR13, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR13out <- scanone(HvDRR13per, method = "hk")
plot(HvDRR13out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR13permout <- scanone(HvDRR13per, method = "hk", n.perm = 4000)
summary(HvDRR13out, perms = HvDRR13permout, alpha = 0.05, pvalues = TRUE)
HvDRR13qtl1 <- makeqtl(HvDRR13per, chr = c(3,6), pos = c(69.4, 72), what = "prob")
HvDRR13outc36 <- addqtl(HvDRR13per, qtl = HvDRR13qtl1, method = "hk")
summary(HvDRR13outc36, perms = HvDRR13permout, alpha = 0.05, pvalues = TRUE)
HvDRR13qtl2 <- makeqtl(HvDRR13per, chr = c(3,6,7), pos = c(69.4, 72, 53.1), what = "prob")
HvDRR13outc367 <- addqtl(HvDRR13per, qtl = HvDRR13qtl2, method = "hk")
summary(HvDRR13outc367, perms = HvDRR13permout, alpha = 0.05, pvalues = TRUE)

HvDRR13rqtl <- refineqtl(HvDRR13per, qtl = HvDRR13qtl2, method = "hk", verbose = FALSE)
summary(HvDRR13rqtl)
summary(fitqtl(HvDRR13per, qtl = HvDRR13rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR13per, 3, 66.7)
effectplot(HvDRR13per, mname1 = "JHI-Hv50k-2016-164723")
summary(HvDRR13per, mname1 = "JHI-Hv50k-2016-164723")
plotPXG(HvDRR13per, "JHI-Hv50k-2016-164723")

find.marker(HvDRR13per, 6, 69)
effectplot(HvDRR13per, mname1 = "JHI-Hv50k-2016-383599")
summary(HvDRR13per, mname1 = "JHI-Hv50k-2016-383599")
plotPXG(HvDRR13per, "JHI-Hv50k-2016-383599")

find.marker(HvDRR13per, 7, 53.1)
effectplot(HvDRR13per, mname1 = "JHI-Hv50k-2016-458766")
summary(HvDRR13per, mname1 = "JHI-Hv50k-2016-458766")
plotPXG(HvDRR13per, "JHI-Hv50k-2016-458766")

summary(fitqtl(HvDRR13per, qtl = HvDRR13rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3))

HvDRR13CI3 <- lodint(HvDRR13rqtl, chr = 3, qtl.index = 1, expandtomarkers = TRUE)
HvDRR13CI6 <- lodint(HvDRR13rqtl, chr = 6, qtl.index = 2, expandtomarkers = TRUE)
HvDRR13CI7 <- lodint(HvDRR13rqtl, chr = 7, qtl.index = 3, expandtomarkers = TRUE)

pt(-4.964, 64, lower.tail = FALSE)
#p<0.01
pt(-4.194, 64, lower.tail = FALSE)
#p<0.01
pt(-4.234, 64, lower.tail = FALSE)
#p<0.01

#HvDRR14####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR14 <- Populations$HvDRR14

HvDRR14$geno$'1'$data<- apply(HvDRR14$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'2'$data<- apply(HvDRR14$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'3'$data<- apply(HvDRR14$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'4'$data<- apply(HvDRR14$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'5'$data<- apply(HvDRR14$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'6'$data<- apply(HvDRR14$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'7'$data<- apply(HvDRR14$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR14$geno$'1'$data<- apply(HvDRR14$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'2'$data<- apply(HvDRR14$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'3'$data<- apply(HvDRR14$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'4'$data<- apply(HvDRR14$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'5'$data<- apply(HvDRR14$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'6'$data<- apply(HvDRR14$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR14$geno$'7'$data<- apply(HvDRR14$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR14)[1] <- "riself"

summary(HvDRR14)
plot(HvDRR14)
HvDRR14 <- subset(HvDRR14,ind=c(1:nrow(HvDRR14$pheno))[!is.na(HvDRR14$pheno)])

HvDRR14per <- calc.genoprob(HvDRR14, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR14out <- scanone(HvDRR14per, method = "hk")
plot(HvDRR14out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR14permout <- scanone(HvDRR14per, method = "hk", n.perm = 4000)
summary(HvDRR14out, perms = HvDRR14permout, alpha = 0.05, pvalues = TRUE)
HvDRR14qtl1 <- makeqtl(HvDRR14per, chr = 3, pos = 178, what = "prob")
HvDRR14outc3 <- addqtl(HvDRR14per, qtl = HvDRR14qtl1, method = "hk")
summary(HvDRR14outc3, perms = HvDRR14permout, alpha = 0.05, pvalues = TRUE)

HvDRR14rqtl <- refineqtl(HvDRR14per, qtl = HvDRR14qtl1, method = "hk", verbose = FALSE)
summary(HvDRR14rqtl)
summary(fitqtl(HvDRR14per, qtl = HvDRR14rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR14per, 3, 178)
effectplot(HvDRR14per, mname1 = "BOPA1_ABC07496-pHv1343-02")
summary(HvDRR14per, mname1 = "BOPA1_ABC07496-pHv1343-02")
plotPXG(HvDRR14per, "BOPA1_ABC07496-pHv1343-02")

summary(fitqtl(HvDRR14per, qtl = HvDRR14rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR14CI3 <- lodint(HvDRR14rqtl, chr = 3, qtl.index = 1, expandtomarkers = TRUE)

pt(4.404, 74, lower.tail = FALSE)
#p<0.01

#HvDRR15####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR15 <- Populations$HvDRR15

HvDRR15$geno$'1'$data<- apply(HvDRR15$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'2'$data<- apply(HvDRR15$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'3'$data<- apply(HvDRR15$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'4'$data<- apply(HvDRR15$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'5'$data<- apply(HvDRR15$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'6'$data<- apply(HvDRR15$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'7'$data<- apply(HvDRR15$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR15$geno$'1'$data<- apply(HvDRR15$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'2'$data<- apply(HvDRR15$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'3'$data<- apply(HvDRR15$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'4'$data<- apply(HvDRR15$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'5'$data<- apply(HvDRR15$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'6'$data<- apply(HvDRR15$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR15$geno$'7'$data<- apply(HvDRR15$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR15)[1] <- "riself"


summary(HvDRR15)
plot(HvDRR15)
HvDRR15 <- subset(HvDRR15,ind=c(1:nrow(HvDRR15$pheno))[!is.na(HvDRR15$pheno)])

HvDRR15per <- calc.genoprob(HvDRR15, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR15out <- scanone(HvDRR15per, method = "hk")
plot(HvDRR15out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR15permout <- scanone(HvDRR15per, method = "hk", n.perm = 4000)
summary(HvDRR15out, perms = HvDRR15permout, alpha = 0.05, pvalues = TRUE)
HvDRR15qtl1 <- makeqtl(HvDRR15per, chr = 3, pos = 98, what = "prob")
HvDRR15outc3 <- addqtl(HvDRR15per, qtl = HvDRR15qtl1, method = "hk")
summary(HvDRR15outc3, perms = HvDRR15permout, alpha = 0.05, pvalues = TRUE)

HvDRR15rqtl <- refineqtl(HvDRR15per, qtl = HvDRR15qtl1, method = "hk", verbose = FALSE)
summary(HvDRR15rqtl)
summary(fitqtl(HvDRR15per, qtl = HvDRR15rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR15per, 3, 98)
effectplot(HvDRR15per, mname1 = "JHI-Hv50k-2016-183463")
summary(HvDRR15per, mname1 = "JHI-Hv50k-2016-183463")
plotPXG(HvDRR15per, "JHI-Hv50k-2016-183463")

summary(fitqtl(HvDRR15per, qtl = HvDRR15rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR15CI3 <- lodint(HvDRR15rqtl, chr = 3, qtl.index = 1, expandtomarkers = TRUE)

pt(-4.242, 45, lower.tail = FALSE)
#p<0.01

#HvDRR16####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR16 <- Populations$HvDRR16

HvDRR16$geno$'1'$data<- apply(HvDRR16$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'2'$data<- apply(HvDRR16$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'3'$data<- apply(HvDRR16$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'4'$data<- apply(HvDRR16$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'5'$data<- apply(HvDRR16$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'6'$data<- apply(HvDRR16$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'7'$data<- apply(HvDRR16$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR16$geno$'1'$data<- apply(HvDRR16$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'2'$data<- apply(HvDRR16$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'3'$data<- apply(HvDRR16$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'4'$data<- apply(HvDRR16$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'5'$data<- apply(HvDRR16$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'6'$data<- apply(HvDRR16$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR16$geno$'7'$data<- apply(HvDRR16$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR16)[1] <- "riself"

summary(HvDRR16)
plot(HvDRR16)
HvDRR16 <- subset(HvDRR16,ind=c(1:nrow(HvDRR16$pheno))[!is.na(HvDRR16$pheno)])

HvDRR16per <- calc.genoprob(HvDRR16, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR16out <- scanone(HvDRR16per, method = "hk")
plot(HvDRR16out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR16permout <- scanone(HvDRR16per, method = "hk", n.perm = 4000)
summary(HvDRR16out, perms = HvDRR16permout, alpha = 0.05, pvalues = TRUE)
HvDRR16qtl1 <- makeqtl(HvDRR16per, chr = 2, pos = 268, what = "prob")
HvDRR16outc2 <- addqtl(HvDRR16per, qtl = HvDRR16qtl1, method = "hk")
summary(HvDRR16outc2, perms = HvDRR16permout, alpha = 0.05, pvalues = TRUE)
HvDRR16qtl2 <- makeqtl(HvDRR16per, c(2, 2), c(268,47), what = "prob")
HvDRR16outc22 <- addqtl(HvDRR16per, qtl = HvDRR16qtl2, method = "hk")
summary(HvDRR16outc22, perms = HvDRR16permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR16out, HvDRR16outc22, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)

HvDRR16rqtl <- refineqtl(HvDRR16per, qtl = HvDRR16qtl2, method = "hk", verbose = FALSE)
summary(HvDRR16rqtl)
summary(fitqtl(HvDRR16per, qtl = HvDRR16rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR16per, 2, 274)
effectplot(HvDRR16per, mname1 = "JHI-Hv50k-2016-131048")
summary(HvDRR16per, mname1 = "JHI-Hv50k-2016-131048")
plotPXG(HvDRR16per, "JHI-Hv50k-2016-131048")

find.marker(HvDRR16per, 2, 47)
effectplot(HvDRR16per, mname1 = "BOPA2_12_30872")
summary(HvDRR16per, mname1 = "BOPA2_12_30872")
plotPXG(HvDRR16per, "BOPA2_12_30872")

summary(fitqtl(HvDRR16per, qtl = HvDRR16rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1+Q2))

HvDRR16CI2 <- lodint(HvDRR16rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR16CI22 <- lodint(HvDRR16rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)

pt(-6.590, 86, lower.tail = FALSE)
#p<0.01
pt(-5.075, 86, lower.tail = FALSE)
#p<0.01

#HvDRR17####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR17 <- Populations$HvDRR17

HvDRR17$geno$'1'$data<- apply(HvDRR17$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'2'$data<- apply(HvDRR17$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'3'$data<- apply(HvDRR17$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'4'$data<- apply(HvDRR17$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'5'$data<- apply(HvDRR17$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'6'$data<- apply(HvDRR17$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'7'$data<- apply(HvDRR17$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR17$geno$'1'$data<- apply(HvDRR17$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'2'$data<- apply(HvDRR17$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'3'$data<- apply(HvDRR17$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'4'$data<- apply(HvDRR17$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'5'$data<- apply(HvDRR17$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'6'$data<- apply(HvDRR17$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR17$geno$'7'$data<- apply(HvDRR17$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR17)[1] <- "riself"

summary(HvDRR17)
plot(HvDRR17)
HvDRR17 <- subset(HvDRR17,ind=c(1:nrow(HvDRR17$pheno))[!is.na(HvDRR17$pheno)])

HvDRR17per <- calc.genoprob(HvDRR17, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR17out <- scanone(HvDRR17per, method = "hk")
plot(HvDRR17out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR17permout <- scanone(HvDRR17per, method = "hk", n.perm=4000)
summary(HvDRR17out, perms = HvDRR17permout, alpha = 0.05, pvalues = TRUE)

#HvDRR18####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR18 <- Populations$HvDRR18

HvDRR18$geno$'1'$data<- apply(HvDRR18$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'2'$data<- apply(HvDRR18$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'3'$data<- apply(HvDRR18$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'4'$data<- apply(HvDRR18$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'5'$data<- apply(HvDRR18$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'6'$data<- apply(HvDRR18$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'7'$data<- apply(HvDRR18$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR18$geno$'1'$data<- apply(HvDRR18$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'2'$data<- apply(HvDRR18$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'3'$data<- apply(HvDRR18$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'4'$data<- apply(HvDRR18$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'5'$data<- apply(HvDRR18$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'6'$data<- apply(HvDRR18$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR18$geno$'7'$data<- apply(HvDRR18$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR18)[1] <- "riself"

summary(HvDRR18)
plot(HvDRR18)
HvDRR18 <- subset(HvDRR18,ind=c(1:nrow(HvDRR18$pheno))[!is.na(HvDRR18$pheno)])

HvDRR18per <- calc.genoprob(HvDRR18, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR18out <- scanone(HvDRR18per, method = "hk")
plot(HvDRR18out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR18permout <- scanone(HvDRR18per, method = "hk", n.perm = 4000)
summary(HvDRR18out, perms = HvDRR18permout, alpha = 0.05, pvalues = TRUE)
HvDRR18qtl1 <- makeqtl(HvDRR18per, chr = 7, pos = 37.1, what = "prob")
HvDRR18outc7 <- addqtl(HvDRR18per, qtl = HvDRR18qtl1, method = "hk")
summary(HvDRR18outc7, perms = HvDRR18permout, alpha = 0.05, pvalues = TRUE)

HvDRR18rqtl <- refineqtl(HvDRR18per, qtl = HvDRR18qtl1, method = "hk", verbose = FALSE)
summary(HvDRR18rqtl)
summary(fitqtl(HvDRR18per, qtl = HvDRR18rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR18per, 7, 37.1)
effectplot(HvDRR18per, mname1 = "JHI-Hv50k-2016-454940")
summary(HvDRR18per, mname1 = "JHI-Hv50k-2016-454940")
plotPXG(HvDRR18per, "JHI-Hv50k-2016-454940")

summary(fitqtl(HvDRR18per, qtl = HvDRR18rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR18CI7 <- lodint(HvDRR18rqtl, chr = 7, qtl.index = 1, expandtomarkers = TRUE)

pt(-4.167, 77, lower.tail = FALSE)
#p<0.01

#HvDRR19####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR19 <- Populations$HvDRR19

HvDRR19$geno$'1'$data<- apply(HvDRR19$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'2'$data<- apply(HvDRR19$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'3'$data<- apply(HvDRR19$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'4'$data<- apply(HvDRR19$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'5'$data<- apply(HvDRR19$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'6'$data<- apply(HvDRR19$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'7'$data<- apply(HvDRR19$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR19$geno$'1'$data<- apply(HvDRR19$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'2'$data<- apply(HvDRR19$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'3'$data<- apply(HvDRR19$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'4'$data<- apply(HvDRR19$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'5'$data<- apply(HvDRR19$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'6'$data<- apply(HvDRR19$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR19$geno$'7'$data<- apply(HvDRR19$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR19)[1] <- "riself"

summary(HvDRR19)
plot(HvDRR19)
HvDRR19 <- subset(HvDRR19,ind=c(1:nrow(HvDRR19$pheno))[!is.na(HvDRR19$pheno)])

HvDRR19per <- calc.genoprob(HvDRR19, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR19out <- scanone(HvDRR19per, method = "hk")
plot(HvDRR19out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR19permout <- scanone(HvDRR19per, method = "hk", n.perm = 4000)
summary(HvDRR19out, perms = HvDRR19permout, alpha = 0.05, pvalues = TRUE)
HvDRR19qtl1 <- makeqtl(HvDRR19per, chr = c(2, 3), pos = c(104,192), what = "prob")
HvDRR19outc23 <- addqtl(HvDRR19per, qtl = HvDRR19qtl1, method = "hk")
summary(HvDRR19outc23, perms = HvDRR19permout, alpha = 0.05, pvalues = TRUE)

HvDRR19rqtl <- refineqtl(HvDRR19per, qtl = HvDRR19qtl1, method = "hk", verbose = FALSE)
summary(HvDRR19rqtl)
summary(fitqtl(HvDRR19per, qtl = HvDRR19rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR19per, 2, 103.6)
effectplot(HvDRR19per, mname1 = "JHI-Hv50k-2016-98990")
summary(HvDRR19per, mname1 = "JHI-Hv50k-2016-98990")
plotPXG(HvDRR19per, "JHI-Hv50k-2016-98990")

find.marker(HvDRR19per, 3, 192.2)
effectplot(HvDRR19per, mname1 = "JHI-Hv50k-2016-203690")
summary(HvDRR19per, mname1 = "JHI-Hv50k-2016-203690")
plotPXG(HvDRR19per, "JHI-Hv50k-2016-203690")

summary(fitqtl(HvDRR19per, qtl = HvDRR19rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1 + Q2))

HvDRR19CI2 <- lodint(HvDRR19rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR19CI3 <- lodint(HvDRR19rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)

pt(-4.200, 80, lower.tail = FALSE)
#p<0.01
pt(-3.764, 80, lower.tail = FALSE)
#p<0.05
#HvDRR20####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR20 <- Populations$HvDRR20

HvDRR20$geno$'1'$data<- apply(HvDRR20$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'2'$data<- apply(HvDRR20$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'3'$data<- apply(HvDRR20$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'4'$data<- apply(HvDRR20$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'5'$data<- apply(HvDRR20$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'6'$data<- apply(HvDRR20$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'7'$data<- apply(HvDRR20$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR20$geno$'1'$data<- apply(HvDRR20$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'2'$data<- apply(HvDRR20$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'3'$data<- apply(HvDRR20$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'4'$data<- apply(HvDRR20$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'5'$data<- apply(HvDRR20$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'6'$data<- apply(HvDRR20$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR20$geno$'7'$data<- apply(HvDRR20$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR20)[1] <- "riself"

summary(HvDRR20)
plot(HvDRR20)
HvDRR20 <- subset(HvDRR20,ind=c(1:nrow(HvDRR20$pheno))[!is.na(HvDRR20$pheno)])

HvDRR20per <- calc.genoprob(HvDRR20, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR20out <- scanone(HvDRR20per, method = "hk")
plot(HvDRR20out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR20permout <- scanone(HvDRR20per, method = "hk", n.perm = 4000)
summary(HvDRR20out, perms = HvDRR20permout, alpha = 0.05, pvalues = TRUE)
HvDRR20qtl1 <- makeqtl(HvDRR20per, chr = 7, pos = 29.6, what = "prob")
HvDRR20outc7 <- addqtl(HvDRR20per, qtl = HvDRR20qtl1, method = "hk")
summary(HvDRR20outc7, perms = HvDRR20permout, alpha = 0.05, pvalues = TRUE)

HvDRR20rqtl <- refineqtl(HvDRR20per, qtl = HvDRR20qtl1, method = "hk", verbose = FALSE)
summary(HvDRR20rqtl)
summary(fitqtl(HvDRR20per, qtl = HvDRR20rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR20per, 7, 29.6)
effectplot(HvDRR20per, mname1 = "JHI-Hv50k-2016-455231")
summary(HvDRR20per, mname1 = "JHI-Hv50k-2016-455231")
plotPXG(HvDRR20per, "JHI-Hv50k-2016-455231")

summary(fitqtl(HvDRR20per, qtl = HvDRR20rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR20CI7 <- lodint(HvDRR20rqtl, chr = 7, qtl.index = 1, expandtomarkers = TRUE)

pt(4.253, 96, lower.tail = FALSE)
#p<0.01

#HvDRR21####

#setwd("/gpfs/project/projects/qggp/Francesco_barley/1stpaper/Rprojects/QTLanalysis/")
Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
#library(qtl)
HvDRR21 <- Populations$HvDRR21

HvDRR21$geno$'1'$data<- apply(HvDRR21$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'2'$data<- apply(HvDRR21$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'3'$data<- apply(HvDRR21$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'4'$data<- apply(HvDRR21$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'5'$data<- apply(HvDRR21$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'6'$data<- apply(HvDRR21$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'7'$data<- apply(HvDRR21$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR21$geno$'1'$data<- apply(HvDRR21$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'2'$data<- apply(HvDRR21$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'3'$data<- apply(HvDRR21$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'4'$data<- apply(HvDRR21$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'5'$data<- apply(HvDRR21$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'6'$data<- apply(HvDRR21$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR21$geno$'7'$data<- apply(HvDRR21$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR21)[1] <- "riself"

summary(HvDRR21)
plot(HvDRR21)
HvDRR21 <- subset(HvDRR21,ind=c(1:nrow(HvDRR21$pheno))[!is.na(HvDRR21$pheno)])

HvDRR21per <- calc.genoprob(HvDRR21, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR21out <- scanone(HvDRR21per, method = "hk")
plot(HvDRR21out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR21permout <- scanone(HvDRR21per, method = "hk", n.perm = 4000)
summary(HvDRR21out, perms = HvDRR21permout, alpha = 0.05, pvalues = TRUE)
HvDRR21qtl1 <- makeqtl(HvDRR21per, chr = 7, pos = 135, what = "prob")
HvDRR21outc7 <- addqtl(HvDRR21per, qtl = HvDRR21qtl1, method = "hk")
summary(HvDRR21outc7, perms = HvDRR21permout, alpha = 0.05, pvalues = TRUE)
HvDRR21qtl2 <- makeqtl(HvDRR21per, chr = c(7,7), pos = c(135,65.5), what = "prob")
HvDRR21outc77 <- addqtl(HvDRR21per, qtl = HvDRR21qtl2, method = "hk")
summary(HvDRR21outc77, perms = HvDRR21permout, alpha = 0.05, pvalues = TRUE)

HvDRR21rqtl <- refineqtl(HvDRR21per, qtl = HvDRR21qtl2, method = "hk", verbose = FALSE)
summary(HvDRR21rqtl)
summary(fitqtl(HvDRR21per, qtl = HvDRR21rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR21per, 7, 135)
effectplot(HvDRR21per, mname1 = "JHI-Hv50k-2016-481738")
summary(HvDRR21per, mname1 = "JHI-Hv50k-2016-481738")
plotPXG(HvDRR21per, "JHI-Hv50k-2016-481738")

find.marker(HvDRR21per, 7, 65.5)
effectplot(HvDRR21per, mname1 = "JHI-Hv50k-2016-459234")
summary(HvDRR21per, mname1 = "JHI-Hv50k-2016-459234")
plotPXG(HvDRR21per, "JHI-Hv50k-2016-459234")

summary(fitqtl(HvDRR21per, qtl = HvDRR21rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1+Q2))

HvDRR21CI7 <- lodint(HvDRR21rqtl, chr = 7, qtl.index = 1, expandtomarkers = TRUE)
HvDRR21CI77 <- lodint(HvDRR21rqtl, chr = 7, qtl.index = 2, expandtomarkers = TRUE)

pt(-6.186, 68, lower.tail = FALSE)
#p<0.01
pt(-3.923, 68, lower.tail = FALSE)
#p<0.01

#HvDRR22####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR22 <- Populations$HvDRR22

HvDRR22$geno$'1'$data<- apply(HvDRR22$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'2'$data<- apply(HvDRR22$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'3'$data<- apply(HvDRR22$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'4'$data<- apply(HvDRR22$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'5'$data<- apply(HvDRR22$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'6'$data<- apply(HvDRR22$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'7'$data<- apply(HvDRR22$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR22$geno$'1'$data<- apply(HvDRR22$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'2'$data<- apply(HvDRR22$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'3'$data<- apply(HvDRR22$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'4'$data<- apply(HvDRR22$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'5'$data<- apply(HvDRR22$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'6'$data<- apply(HvDRR22$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR22$geno$'7'$data<- apply(HvDRR22$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR22)[1] <- "riself"

summary(HvDRR22)
plot(HvDRR22)
HvDRR22 <- subset(HvDRR22,ind=c(1:nrow(HvDRR22$pheno))[!is.na(HvDRR22$pheno)])

HvDRR22per <- calc.genoprob(HvDRR22, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR22out <- scanone(HvDRR22per, method = "hk")
plot(HvDRR22out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR22permout <- scanone(HvDRR22per, method = "hk", n.perm=4000)
summary(HvDRR22out, perms = HvDRR22permout, alpha = 0.05, pvalues = TRUE)
HvDRR22qtl1 <- makeqtl(HvDRR22per, chr = 7, pos = 119, what = "prob")
HvDRR22outc7 <- addqtl(HvDRR22per, qtl = HvDRR22qtl1, method = "hk")
summary(HvDRR22outc7, perms = HvDRR22permout, alpha = 0.05, pvalues = TRUE)
HvDRR22qtl2 <- makeqtl(HvDRR22per, chr = c(7,1,2,7), pos = c(119,217.6,84.7,64), what = "prob")
HvDRR22outc7127 <- addqtl(HvDRR22per, qtl = HvDRR22qtl2, method = "hk")
summary(HvDRR22outc7127, perms = HvDRR22permout, alpha = 0.05, pvalues = TRUE)
HvDRR22qtl3 <- makeqtl(HvDRR22per, chr = c(7,1,2,7,4), pos = c(119,217.6,84.7,64,65), what = "prob")
HvDRR22outc71274 <- addqtl(HvDRR22per, qtl = HvDRR22qtl3, method = "hk")
summary(HvDRR22outc71274, perms = HvDRR22permout, alpha = 0.05, pvalues = TRUE)
HvDRR22qtl4 <- makeqtl(HvDRR22per, chr = c(7,1,2,7,4,6), pos = c(119,217.6,84.7,64,65,57.5), what = "prob")
HvDRR22outc712746 <- addqtl(HvDRR22per, qtl = HvDRR22qtl4, method = "hk")
summary(HvDRR22outc712746, perms = HvDRR22permout, alpha = 0.05, pvalues = TRUE)

HvDRR22rqtl <- refineqtl(HvDRR22per, qtl = HvDRR22qtl4, method = "hk", verbose = FALSE)
summary(HvDRR22rqtl)
summary(fitqtl(HvDRR22per, qtl = HvDRR22rqtl, method = "hk"), pvalues = FALSE)


find.marker(HvDRR22per, 7, 119.3)
effectplot(HvDRR22per, mname1 = "JHI-Hv50k-2016-479764")
summary(HvDRR22per, mname1 = "JHI-Hv50k-2016-479764")
plotPXG(HvDRR22per, "JHI-Hv50k-2016-479764")

find.marker(HvDRR22per, 1, 214)
effectplot(HvDRR22per, mname1 = "SCRI_RS_147318")
summary(HvDRR22per, mname1 = "SCRI_RS_147318")
plotPXG(HvDRR22per, "SCRI_RS_147318")

find.marker(HvDRR22per, 2, 85.2)
effectplot(HvDRR22per, mname1 = "JHI-Hv50k-2016-93257")
summary(HvDRR22per, mname1 = "JHI-Hv50k-2016-93257")
plotPXG(HvDRR22per, "JHI-Hv50k-2016-93257")

find.marker(HvDRR22per, 7, 63)
effectplot(HvDRR22per, mname1 = "JHI-Hv50k-2016-460104")
summary(HvDRR22per, mname1 = "JHI-Hv50k-2016-460104")
plotPXG(HvDRR22per, "JHI-Hv50k-2016-460104")

find.marker(HvDRR22per, 4, 64.4)
effectplot(HvDRR22per, mname1 = "JHI-Hv50k-2016-242645")
summary(HvDRR22per, mname1 = "JHI-Hv50k-2016-242645")
plotPXG(HvDRR22per, "JHI-Hv50k-2016-242645")

find.marker(HvDRR22per, 6, 57.5)
effectplot(HvDRR22per, mname1 = "SCRI_RS_118381")
summary(HvDRR22per, mname1 = "SCRI_RS_118381")
plotPXG(HvDRR22per, "SCRI_RS_118381")

summary(fitqtl(HvDRR22per, qtl = HvDRR22rqtl, method = "hk", get.ests = TRUE, formula = Flowering ~Q1 + Q2 + Q3 + Q4 + Q5 + Q6))

HvDRR22CI7 <- lodint(HvDRR22rqtl, chr = 7, qtl.index = 1, expandtomarkers = TRUE)
HvDRR22CI1 <- lodint(HvDRR22rqtl, chr = 1, qtl.index = 2, expandtomarkers = TRUE)
HvDRR22CI2 <- lodint(HvDRR22rqtl, chr = 2, qtl.index = 3, expandtomarkers = TRUE)
HvDRR22CI77 <- lodint(HvDRR22rqtl, chr = 7, qtl.index = 4, expandtomarkers = TRUE)
HvDRR22CI4 <- lodint(HvDRR22rqtl, chr = 4, qtl.index = 5, expandtomarkers = TRUE)
HvDRR22CI6 <- lodint(HvDRR22rqtl, chr = 6, qtl.index = 6, expandtomarkers = TRUE)

pt(-12.544, 83, lower.tail = FALSE)
#p<0.01
pt(-3.311, 83, lower.tail = FALSE)
#p<0.01
pt(-5.389, 83, lower.tail = FALSE)
#p<0.01
pt(-6.261, 83, lower.tail = FALSE)
#p<0.01
pt(5.661, 83, lower.tail = FALSE)
#p<0.01
pt(4.543, 83, lower.tail = FALSE)
#p<0.01

#HvDRR23####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR23 <- Populations$HvDRR23

HvDRR23$geno$'1'$data<- apply(HvDRR23$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'2'$data<- apply(HvDRR23$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'3'$data<- apply(HvDRR23$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'4'$data<- apply(HvDRR23$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'5'$data<- apply(HvDRR23$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'6'$data<- apply(HvDRR23$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'7'$data<- apply(HvDRR23$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR23$geno$'1'$data<- apply(HvDRR23$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'2'$data<- apply(HvDRR23$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'3'$data<- apply(HvDRR23$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'4'$data<- apply(HvDRR23$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'5'$data<- apply(HvDRR23$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'6'$data<- apply(HvDRR23$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR23$geno$'7'$data<- apply(HvDRR23$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR23)[1] <- "riself"

summary(HvDRR23)
plot(HvDRR23)
HvDRR23 <- subset(HvDRR23,ind=c(1:nrow(HvDRR23$pheno))[!is.na(HvDRR23$pheno)])

HvDRR23per <- calc.genoprob(HvDRR23, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR23out <- scanone(HvDRR23per, method = "hk")
plot(HvDRR23out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR23permout <- scanone(HvDRR23per, method = "hk", n.perm = 4000)
summary(HvDRR23out, perms = HvDRR23permout, alpha = 0.05, pvalues = TRUE)
HvDRR23qtl1 <- makeqtl(HvDRR23per, chr = 7, pos = 116, what = "prob")
HvDRR23outc7 <- addqtl(HvDRR23per, qtl = HvDRR23qtl1, method = "hk")
summary(HvDRR23outc7, perms = HvDRR23permout, alpha = 0.05, pvalues = TRUE)

HvDRR23rqtl <- refineqtl(HvDRR23per, qtl = HvDRR23qtl1, method = "hk", verbose = FALSE)
summary(HvDRR23rqtl)
summary(fitqtl(HvDRR23per, qtl = HvDRR23rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR23per, 7, 116)
effectplot(HvDRR23per, mname1 = "JHI-Hv50k-2016-479764")
summary(HvDRR23per, mname1 = "JHI-Hv50k-2016-479764")
plotPXG(HvDRR23per, "JHI-Hv50k-2016-479764")

summary(fitqtl(HvDRR23per, qtl = HvDRR23rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1))

HvDRR23CI7 <- lodint(HvDRR23rqtl, chr = 7, qtl.index = 1, expandtomarkers = TRUE)

pt(-6.343, 84, lower.tail = FALSE)
#p<0.01
#HvDRR24####

#setwd("/gpfs/project/projects/qggp/Francesco_barley/1stpaper/Rprojects/QTLanalysis/")
Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
#library(qtl)
HvDRR24 <- Populations$HvDRR24

HvDRR24$geno$'1'$data<- apply(HvDRR24$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'2'$data<- apply(HvDRR24$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'3'$data<- apply(HvDRR24$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'4'$data<- apply(HvDRR24$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'5'$data<- apply(HvDRR24$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'6'$data<- apply(HvDRR24$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'7'$data<- apply(HvDRR24$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR24$geno$'1'$data<- apply(HvDRR24$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'2'$data<- apply(HvDRR24$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'3'$data<- apply(HvDRR24$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'4'$data<- apply(HvDRR24$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'5'$data<- apply(HvDRR24$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'6'$data<- apply(HvDRR24$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR24$geno$'7'$data<- apply(HvDRR24$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR24)[1] <- "riself"

summary(HvDRR24)
plot(HvDRR24)
HvDRR24 <- subset(HvDRR24,ind=c(1:nrow(HvDRR24$pheno))[!is.na(HvDRR24$pheno)])

HvDRR24per <- calc.genoprob(HvDRR24, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR24out <- scanone(HvDRR24per, method = "hk")
plot(HvDRR24out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR24permout <- scanone(HvDRR24per, method = "hk", n.perm=4000)
summary(HvDRR24out, perms = HvDRR24permout, alpha = 0.05, pvalues = TRUE)
HvDRR24qtl1 <- makeqtl(HvDRR24per, chr = 3, pos = 20, what = "prob")
HvDRR24outc3 <- addqtl(HvDRR24per, qtl = HvDRR24qtl1, method = "hk")
summary(HvDRR24outc3, perms = HvDRR24permout, alpha = 0.05, pvalues = TRUE)
HvDRR24qtl2 <- makeqtl(HvDRR24per, chr = c(3,1), pos = c(20,106), what = "prob")
HvDRR24outc31 <- addqtl(HvDRR24per, qtl = HvDRR24qtl2, method = "hk")
summary(HvDRR24outc31, perms = HvDRR24permout, alpha = 0.05, pvalues = TRUE)

HvDRR24rqtl <- refineqtl(HvDRR24per, qtl = HvDRR24qtl2, method = "hk", verbose = FALSE)
summary(HvDRR24rqtl)
summary(fitqtl(HvDRR24per, qtl = HvDRR24rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR24per, 3, 20)
effectplot(HvDRR24per, mname1 = "JHI-Hv50k-2016-154007")
summary(HvDRR24per, mname1 = "JHI-Hv50k-2016-154007")
plotPXG(HvDRR24per, "JHI-Hv50k-2016-154007")

find.marker(HvDRR24per, 1, 106)
effectplot(HvDRR24per, mname1 = "JHI-Hv50k-2016-34594")
summary(HvDRR24per, mname1 = "JHI-Hv50k-2016-34594")
plotPXG(HvDRR24per, "JHI-Hv50k-2016-34594")

summary(fitqtl(HvDRR24per, qtl = HvDRR24rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2))

HvDRR24CI3 <- lodint(HvDRR24rqtl, chr = 3, qtl.index = 1, expandtomarkers = TRUE)
HvDRR24CI1 <- lodint(HvDRR24rqtl, chr = 1, qtl.index = 2, expandtomarkers = TRUE)

pt(-5.210, 61, lower.tail = FALSE)
#p<0.01
pt(-4.129, 61, lower.tail = FALSE)
#p<0.01

#HvDRR25####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR25 <- Populations$HvDRR25

HvDRR25$geno$'1'$data<- apply(HvDRR25$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'2'$data<- apply(HvDRR25$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'3'$data<- apply(HvDRR25$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'4'$data<- apply(HvDRR25$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'5'$data<- apply(HvDRR25$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'6'$data<- apply(HvDRR25$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR25$geno$'7'$data<- apply(HvDRR25$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
#to eliminate heterozygous
HvDRR25$geno$'1'$data<- apply(HvDRR25$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'2'$data<- apply(HvDRR25$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'3'$data<- apply(HvDRR25$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'4'$data<- apply(HvDRR25$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'5'$data<- apply(HvDRR25$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'6'$data<- apply(HvDRR25$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR25$geno$'7'$data<- apply(HvDRR25$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR25)[1] <- "riself"


summary(HvDRR25)
plot(HvDRR25)
HvDRR25 <- subset(HvDRR25,ind=c(1:nrow(HvDRR25$pheno))[!is.na(HvDRR25$pheno)])

testmarkername <- as.data.frame(HvDRR25$geno$`3`$map)
which(row.names(testmarkername) == 'BOPA1_5260-462')
marker2drop <- markernames(HvDRR25, chr = 3)[145]
HvDRR25 <- drop.markers(HvDRR25, marker2drop)

HvDRR25per <- calc.genoprob(HvDRR25, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR25out <- scanone(HvDRR25per, method = "hk")
plot(HvDRR25out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR25permout <- scanone(HvDRR25per, method = "hk", n.perm=4000)
summary(HvDRR25out, perms = HvDRR25permout, alpha = 0.05, pvalues = TRUE)
HvDRR25qtl1 <- makeqtl(HvDRR25per, chr = c(2,3), pos = c(33,151), what = "prob")
HvDRR25outc2 <- addqtl(HvDRR25per, qtl = HvDRR25qtl1, method = "hk")
summary(HvDRR25outc2, perms = HvDRR25permout, alpha = 0.05, pvalues = TRUE)

HvDRR25rqtl <- refineqtl(HvDRR25per, qtl = HvDRR25qtl1, method = "hk", verbose = FALSE)
summary(HvDRR25rqtl)
summary(fitqtl(HvDRR25per, qtl = HvDRR25rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR25per, 2,48.8)
effectplot(HvDRR25per, mname1 = "BOPA1_2646-1277")
summary(HvDRR25per, mname1 = "BOPA1_2646-1277")
plotPXG(HvDRR25per, "BOPA1_2646-1277")

find.marker(HvDRR25per, 3,152)
effectplot(HvDRR25per, mname1 = "JHI-Hv50k-2016-205406")
summary(HvDRR25per, mname1 = "JHI-Hv50k-2016-205406")
plotPXG(HvDRR25per, "JHI-Hv50k-2016-205406")

summary(fitqtl(HvDRR25per, qtl = HvDRR25rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2))

HvDRR25CI2 <- lodint(HvDRR25rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR25CI3 <- lodint(HvDRR25rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)

pt(-5.380, 71, lower.tail = FALSE)
pt(4.821, 71, lower.tail = FALSE)

#HvDRR26####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR26 <- Populations$HvDRR26

HvDRR26$geno$'1'$data<- apply(HvDRR26$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'2'$data<- apply(HvDRR26$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'3'$data<- apply(HvDRR26$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'4'$data<- apply(HvDRR26$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'5'$data<- apply(HvDRR26$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'6'$data<- apply(HvDRR26$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'7'$data<- apply(HvDRR26$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR26$geno$'1'$data<- apply(HvDRR26$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'2'$data<- apply(HvDRR26$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'3'$data<- apply(HvDRR26$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'4'$data<- apply(HvDRR26$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'5'$data<- apply(HvDRR26$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'6'$data<- apply(HvDRR26$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR26$geno$'7'$data<- apply(HvDRR26$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR26)[1] <- "riself"

summary(HvDRR26)
plot(HvDRR26)
HvDRR26 <- subset(HvDRR26,ind=c(1:nrow(HvDRR26$pheno))[!is.na(HvDRR26$pheno)])


HvDRR26per <- calc.genoprob(HvDRR26, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR26out <- scanone(HvDRR26per, method = "hk")
plot(HvDRR26out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR26permout <- scanone(HvDRR26per, method = "hk", n.perm=4000)
summary(HvDRR26out, perms = HvDRR26permout, alpha = 0.05, pvalues = TRUE)
HvDRR26qtl1 <- makeqtl(HvDRR26per, chr = c(2,3,4), pos = c(199,86.6,59), what = "prob")
HvDRR26outc234 <- addqtl(HvDRR26per, qtl = HvDRR26qtl1, method = "hk")
summary(HvDRR26outc234, perms = HvDRR26permout, alpha = 0.05, pvalues = TRUE)
HvDRR26qtl2 <- makeqtl(HvDRR26per, chr = c(2,3,4,2), pos = c(199,86.6,59,161), what = "prob")
HvDRR26outc2342 <- addqtl(HvDRR26per, qtl = HvDRR26qtl2, method = "hk")
summary(HvDRR26outc2342, perms = HvDRR26permout, alpha = 0.05, pvalues = TRUE)

HvDRR26rqtl <- refineqtl(HvDRR26per, qtl = HvDRR26qtl2, method = "hk", verbose = FALSE)
summary(HvDRR26rqtl)
summary(fitqtl(HvDRR26per, qtl = HvDRR26rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR26per, 2, 218)
effectplot(HvDRR26per, mname1 = "BOPA2_12_30396")
summary(HvDRR26per, mname1 = "BOPA2_12_30396")
plotPXG(HvDRR26per, "BOPA2_12_30396")

find.marker(HvDRR26per, 3, 87.4)
effectplot(HvDRR26per, mname1 = "JHI-Hv50k-2016-186194")
summary(HvDRR26per, mname1 = "JHI-Hv50k-2016-186194")
plotPXG(HvDRR26per, "JHI-Hv50k-2016-186194")

find.marker(HvDRR26per, 4, 67.4)
effectplot(HvDRR26per, mname1 = "JHI-Hv50k-2016-232799")
summary(HvDRR26per, mname1 = "JHI-Hv50k-2016-232799")
plotPXG(HvDRR26per, "JHI-Hv50k-2016-232799")

find.marker(HvDRR26per, 2, 161.3)
effectplot(HvDRR26per, mname1 = "JHI-Hv50k-2016-109198")
summary(HvDRR26per, mname1 = "JHI-Hv50k-2016-109198")
plotPXG(HvDRR26per, "JHI-Hv50k-2016-109198")

summary(fitqtl(HvDRR26per, qtl = HvDRR26rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3 + Q4))

HvDRR26CI2 <- lodint(HvDRR26rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR26CI3 <- lodint(HvDRR26rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)
HvDRR26CI4 <- lodint(HvDRR26rqtl, chr = 4, qtl.index = 3, expandtomarkers = TRUE)
HvDRR26CI22 <- lodint(HvDRR26rqtl, chr = 2, qtl.index = 4, expandtomarkers = TRUE)

pt(-5.577, 49, lower.tail = FALSE)
#p<0.01
pt(5.055, 49, lower.tail = FALSE)
#p<0.01
pt(5.332, 49, lower.tail = FALSE)
#p<0.01
pt(3.375, 49, lower.tail = FALSE)
#p<0.01

#HvDRR27####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR27 <- Populations$HvDRR27

HvDRR27$geno$'1'$data<- apply(HvDRR27$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'2'$data<- apply(HvDRR27$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'3'$data<- apply(HvDRR27$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'4'$data<- apply(HvDRR27$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'5'$data<- apply(HvDRR27$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'6'$data<- apply(HvDRR27$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'7'$data<- apply(HvDRR27$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR27$geno$'1'$data<- apply(HvDRR27$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'2'$data<- apply(HvDRR27$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'3'$data<- apply(HvDRR27$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'4'$data<- apply(HvDRR27$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'5'$data<- apply(HvDRR27$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'6'$data<- apply(HvDRR27$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR27$geno$'7'$data<- apply(HvDRR27$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR27)[1] <- "riself"

summary(HvDRR27)
plot(HvDRR27)
HvDRR27 <- subset(HvDRR27,ind=c(1:nrow(HvDRR27$pheno))[!is.na(HvDRR27$pheno)])

HvDRR27per <- calc.genoprob(HvDRR27, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR27out <- scanone(HvDRR27per, method = "hk")
plot(HvDRR27out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR27permout <- scanone(HvDRR27per, method = "hk", n.perm=4000)
summary(HvDRR27out, perms = HvDRR27permout, alpha = 0.05, pvalues = TRUE)
HvDRR27qtl1 <- makeqtl(HvDRR27per, chr = 2, pos = 33, what = "prob")
HvDRR27outc2 <- addqtl(HvDRR27per, qtl = HvDRR27qtl1, method = "hk")
summary(HvDRR27outc2, perms = HvDRR27permout, alpha = 0.05, pvalues = TRUE)

HvDRR27rqtl <- refineqtl(HvDRR27per, qtl = HvDRR27qtl1, method = "hk", verbose = FALSE)
summary(HvDRR27rqtl)
summary(fitqtl(HvDRR27per, qtl = HvDRR27rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR27per, 2, 33)
effectplot(HvDRR27per, mname1 = "JHI-Hv50k-2016-71249")
summary(HvDRR27per, mname1 = "JHI-Hv50k-2016-71249")
plotPXG(HvDRR27per, "JHI-Hv50k-2016-71249")

summary(fitqtl(HvDRR27per, qtl = HvDRR27rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR27CI2 <- lodint(HvDRR27rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(4.703, 91, lower.tail = FALSE)
#p<0.01

#HvDRR28####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR28 <- Populations$HvDRR28

HvDRR28$geno$'1'$data<- apply(HvDRR28$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'2'$data<- apply(HvDRR28$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'3'$data<- apply(HvDRR28$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'4'$data<- apply(HvDRR28$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'5'$data<- apply(HvDRR28$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'6'$data<- apply(HvDRR28$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'7'$data<- apply(HvDRR28$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR28$geno$'1'$data<- apply(HvDRR28$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'2'$data<- apply(HvDRR28$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'3'$data<- apply(HvDRR28$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'4'$data<- apply(HvDRR28$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'5'$data<- apply(HvDRR28$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'6'$data<- apply(HvDRR28$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR28$geno$'7'$data<- apply(HvDRR28$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR28)[1] <- "riself"

summary(HvDRR28)
plot(HvDRR28)
HvDRR28 <- subset(HvDRR28,ind=c(1:nrow(HvDRR28$pheno))[!is.na(HvDRR28$pheno)])

HvDRR28per <- calc.genoprob(HvDRR28, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR28out <- scanone(HvDRR28per, method = "hk")
plot(HvDRR28out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR28permout <- scanone(HvDRR28per, method = "hk", n.perm=4000)
summary(HvDRR28out, perms = HvDRR28permout, alpha = 0.05, pvalues = TRUE)

HvDRR28qtl1 <- makeqtl(HvDRR28per, chr = c(2,6), pos = c(47.2,89.3), what = "prob")
HvDRR28outc26 <- addqtl(HvDRR28per, qtl = HvDRR28qtl1, method = "hk")
summary(HvDRR28outc26, perms = HvDRR28permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR28out, HvDRR28outc26, chr = -2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR28out, HvDRR28outc26, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)

HvDRR28rqtl <- refineqtl(HvDRR28per, qtl = HvDRR28qtl1, method = "hk", verbose = FALSE)
summary(HvDRR28rqtl)
summary(fitqtl(HvDRR28per, qtl = HvDRR28rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR28per, 2, 45.1)
effectplot(HvDRR28per, mname1 = "JHI-Hv50k-2016-73692")
summary(HvDRR28per, mname1 = "JHI-Hv50k-2016-73692")
plotPXG(HvDRR28per, "JHI-Hv50k-2016-73692")

find.marker(HvDRR28per, 6, 89.3)
effectplot(HvDRR28per, mname1 = "JHI-Hv50k-2016-408338")
summary(HvDRR28per, mname1 = "JHI-Hv50k-2016-408338")
plotPXG(HvDRR28per, "JHI-Hv50k-2016-408338")

summary(fitqtl(HvDRR28per, qtl = HvDRR28rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2))

HvDRR28CI2 <- lodint(HvDRR28rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR28CI6 <- lodint(HvDRR28rqtl, chr = 6, qtl.index = 2, expandtomarkers = TRUE)

pt(8.033, 78, lower.tail = FALSE)
#p<0.01
pt(-3.954, 78, lower.tail = FALSE)
#p<0.01


#HvDRR29####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR29 <- Populations$HvDRR29

HvDRR29$geno$'1'$data<- apply(HvDRR29$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'2'$data<- apply(HvDRR29$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'3'$data<- apply(HvDRR29$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'4'$data<- apply(HvDRR29$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'5'$data<- apply(HvDRR29$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'6'$data<- apply(HvDRR29$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'7'$data<- apply(HvDRR29$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR29$geno$'1'$data<- apply(HvDRR29$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'2'$data<- apply(HvDRR29$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'3'$data<- apply(HvDRR29$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'4'$data<- apply(HvDRR29$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'5'$data<- apply(HvDRR29$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'6'$data<- apply(HvDRR29$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR29$geno$'7'$data<- apply(HvDRR29$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR29)[1] <- "riself"

summary(HvDRR29)
plot(HvDRR29)
HvDRR29 <- subset(HvDRR29,ind=c(1:nrow(HvDRR29$pheno))[!is.na(HvDRR29$pheno)])

HvDRR29per <- calc.genoprob(HvDRR29, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR29out <- scanone(HvDRR29per, method = "hk")
plot(HvDRR29out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR29permout <- scanone(HvDRR29per, method = "hk", n.perm=4000)
summary(HvDRR29out, perms = HvDRR29permout, alpha = 0.05, pvalues = TRUE)
HvDRR29qtl1 <- makeqtl(HvDRR29per, chr = 2, pos = 146, what = "prob")
HvDRR29outc2 <- addqtl(HvDRR29per, qtl = HvDRR29qtl1, method = "hk")
summary(HvDRR29outc2, perms = HvDRR29permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR29out, HvDRR29outc2, chr = -2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR29out, HvDRR29outc2, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)
HvDRR29qtl2 <- makeqtl(HvDRR29per, c(2,2,6), c(146,43,67), what = "prob")
HvDRR29outc226 <- addqtl(HvDRR29per, qtl = HvDRR29qtl2, method = "hk")
summary(HvDRR29outc226, perms = HvDRR29permout, alpha = 0.05, pvalues = TRUE)

HvDRR29rqtl <- refineqtl(HvDRR29per, qtl = HvDRR29qtl2, method = "hk", verbose = FALSE)
summary(HvDRR29rqtl)
summary(fitqtl(HvDRR29per, qtl = HvDRR29rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR29per, 2, 146.2)
effectplot(HvDRR29per, mname1 = "BOPA2_12_30897")
summary(HvDRR29per, mname1 = "BOPA2_12_30897")
plotPXG(HvDRR29per, "BOPA2_12_30897")

find.marker(HvDRR29per, 2, 43)
effectplot(HvDRR29per, mname1 = "JHI-Hv50k-2016-73780")
summary(HvDRR29per, mname1 = "JHI-Hv50k-2016-73780")
plotPXG(HvDRR29per, "JHI-Hv50k-2016-73780")

find.marker(HvDRR29per, 6, 74)
effectplot(HvDRR29per, mname1 = "JHI-Hv50k-2016-395685")
summary(HvDRR29per, mname1 = "JHI-Hv50k-2016-395685")
plotPXG(HvDRR29per, "JHI-Hv50k-2016-395685")

summary(fitqtl(HvDRR29per, qtl = HvDRR29rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3))

HvDRR29CI2 <- lodint(HvDRR29rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR29CI22 <- lodint(HvDRR29rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)
HvDRR29CI6 <- lodint(HvDRR29rqtl, chr = 6, qtl.index = 3, expandtomarkers = TRUE)

pt(-6.915, 107, lower.tail = FALSE)
#p<0.01
pt(-5.458, 107, lower.tail = FALSE)
#p<0.01
pt(4.275, 107, lower.tail = FALSE)
#p<0.01

#HvDRR30####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR30 <- Populations$HvDRR30

HvDRR30$geno$'1'$data<- apply(HvDRR30$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'2'$data<- apply(HvDRR30$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'3'$data<- apply(HvDRR30$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'4'$data<- apply(HvDRR30$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'5'$data<- apply(HvDRR30$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'6'$data<- apply(HvDRR30$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'7'$data<- apply(HvDRR30$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR30$geno$'1'$data<- apply(HvDRR30$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'2'$data<- apply(HvDRR30$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'3'$data<- apply(HvDRR30$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'4'$data<- apply(HvDRR30$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'5'$data<- apply(HvDRR30$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'6'$data<- apply(HvDRR30$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR30$geno$'7'$data<- apply(HvDRR30$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR30)[1] <- "riself"

summary(HvDRR30)
plot(HvDRR30)
HvDRR30 <- subset(HvDRR30,ind=c(1:nrow(HvDRR30$pheno))[!is.na(HvDRR30$pheno)])

HvDRR30per <- calc.genoprob(HvDRR30, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR30out <- scanone(HvDRR30per, method = "hk")
plot(HvDRR30out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR30permout <- scanone(HvDRR30per, method = "hk", n.perm=4000)
summary(HvDRR30out, perms = HvDRR30permout, alpha = 0.05, pvalues = TRUE)
HvDRR30qtl1 <- makeqtl(HvDRR30per, chr = c(2,4), pos = c(25.7,49), what = "prob")
HvDRR30outc24 <- addqtl(HvDRR30per, qtl = HvDRR30qtl1, method = "hk")
summary(HvDRR30outc24, perms = HvDRR30permout, alpha = 0.05, pvalues = TRUE)

HvDRR30rqtl <- refineqtl(HvDRR30per, qtl = HvDRR30qtl1, method = "hk", verbose = FALSE)
summary(HvDRR30rqtl)
summary(fitqtl(HvDRR30per, qtl = HvDRR30rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR30per, 2, 20)
effectplot(HvDRR30per, mname1 = "JHI-Hv50k-2016-67694")
summary(HvDRR30per, mname1 = "JHI-Hv50k-2016-67694")
plotPXG(HvDRR30per, "JHI-Hv50k-2016-67694")

find.marker(HvDRR30per, 4, 49)
effectplot(HvDRR30per, mname1 = "SCRI_RS_98443")
summary(HvDRR30per, mname1 = "SCRI_RS_98443")
plotPXG(HvDRR30per, "SCRI_RS_98443")

summary(fitqtl(HvDRR30per, qtl = HvDRR30rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2))

HvDRR30CI2 <- lodint(HvDRR30rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR30CI4 <- lodint(HvDRR30rqtl, chr = 4, qtl.index = 2, expandtomarkers = TRUE)

pt(4.571, 120, lower.tail = FALSE)
#p<0.01
pt(4.888, 120, lower.tail = FALSE)
#p<0.01

#HvDRR31####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR31 <- Populations$HvDRR31

HvDRR31$geno$'1'$data<- apply(HvDRR31$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'2'$data<- apply(HvDRR31$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'3'$data<- apply(HvDRR31$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'4'$data<- apply(HvDRR31$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'5'$data<- apply(HvDRR31$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'6'$data<- apply(HvDRR31$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'7'$data<- apply(HvDRR31$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR31$geno$'1'$data<- apply(HvDRR31$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'2'$data<- apply(HvDRR31$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'3'$data<- apply(HvDRR31$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'4'$data<- apply(HvDRR31$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'5'$data<- apply(HvDRR31$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'6'$data<- apply(HvDRR31$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR31$geno$'7'$data<- apply(HvDRR31$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR31)[1] <- "riself"

summary(HvDRR31)
plot(HvDRR31)
HvDRR31 <- subset(HvDRR31,ind=c(1:nrow(HvDRR31$pheno))[!is.na(HvDRR31$pheno)])

HvDRR31per <- calc.genoprob(HvDRR31, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR31out <- scanone(HvDRR31per, method = "hk")
plot(HvDRR31out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR31permout <- scanone(HvDRR31per, method = "hk", n.perm=4000)
summary(HvDRR31out, perms = HvDRR31permout, alpha = 0.05, pvalues = TRUE)
HvDRR31qtl1 <- makeqtl(HvDRR31per, chr = 1, pos = 268, what = "prob")
HvDRR31outc1 <- addqtl(HvDRR31per, qtl = HvDRR31qtl1, method = "hk")
summary(HvDRR31outc1, perms = HvDRR31permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR31out, HvDRR31outc1, chr = -1, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR31out, HvDRR31outc1, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)
HvDRR31qtl2 <- makeqtl(HvDRR31per, chr = c(1,5), pos = c(268,202), what = "prob")
HvDRR31outc15 <- addqtl(HvDRR31per, qtl = HvDRR31qtl2, method = "hk")
summary(HvDRR31outc15, perms = HvDRR31permout, alpha = 0.05, pvalues = TRUE)

HvDRR31rqtl <- refineqtl(HvDRR31per, qtl = HvDRR31qtl2, method = "hk", verbose = FALSE)
summary(HvDRR31rqtl)
summary(fitqtl(HvDRR31per, qtl = HvDRR31rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR31per, 1, 269)
effectplot(HvDRR31per, mname1 = "JHI-Hv50k-2016-58045")
summary(HvDRR31per, mname1 = "JHI-Hv50k-2016-58045")
plotPXG(HvDRR31per, "JHI-Hv50k-2016-58045")

find.marker(HvDRR31per, 5, 202)
effectplot(HvDRR31per, mname1 = "JHI-Hv50k-2016-335543")
summary(HvDRR31per, mname1 = "JHI-Hv50k-2016-335543")
plotPXG(HvDRR31per, "JHI-Hv50k-2016-335543")

summary(fitqtl(HvDRR31per, qtl = HvDRR31rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1 + Q2))

HvDRR31CI1 <- lodint(HvDRR31rqtl, chr = 1, qtl.index = 1, expandtomarkers = TRUE)
HvDRR31CI5 <- lodint(HvDRR31rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)

pt(5.303, 126, lower.tail = FALSE)
pt(4.047, 126, lower.tail = FALSE)

#HvDRR32####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR32 <- Populations$HvDRR32

HvDRR32$geno$'1'$data<- apply(HvDRR32$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'2'$data<- apply(HvDRR32$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'3'$data<- apply(HvDRR32$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'4'$data<- apply(HvDRR32$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'5'$data<- apply(HvDRR32$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'6'$data<- apply(HvDRR32$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'7'$data<- apply(HvDRR32$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR32$geno$'1'$data<- apply(HvDRR32$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'2'$data<- apply(HvDRR32$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'3'$data<- apply(HvDRR32$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'4'$data<- apply(HvDRR32$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'5'$data<- apply(HvDRR32$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'6'$data<- apply(HvDRR32$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR32$geno$'7'$data<- apply(HvDRR32$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR32)[1] <- "riself"

summary(HvDRR32)
plot(HvDRR32)
HvDRR32 <- subset(HvDRR32,ind=c(1:nrow(HvDRR32$pheno))[!is.na(HvDRR32$pheno)])
HvDRR32per <- calc.genoprob(HvDRR32, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR32out <- scanone(HvDRR32per, method = "hk")
plot(HvDRR32out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR32permout <- scanone(HvDRR32per, method = "hk", n.perm=4000)
summary(HvDRR32out, perms = HvDRR32permout, alpha = 0.05, pvalues = TRUE)

#HvDRR33####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR33 <- Populations$HvDRR33

HvDRR33$geno$'1'$data<- apply(HvDRR33$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'2'$data<- apply(HvDRR33$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'3'$data<- apply(HvDRR33$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'4'$data<- apply(HvDRR33$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'5'$data<- apply(HvDRR33$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'6'$data<- apply(HvDRR33$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR33$geno$'7'$data<- apply(HvDRR33$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

HvDRR33$geno$'1'$data<- apply(HvDRR33$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'2'$data<- apply(HvDRR33$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'3'$data<- apply(HvDRR33$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'4'$data<- apply(HvDRR33$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'5'$data<- apply(HvDRR33$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'6'$data<- apply(HvDRR33$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR33$geno$'7'$data<- apply(HvDRR33$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR33)[1] <- "riself"

summary(HvDRR33)
plot(HvDRR33)
HvDRR33 <- subset(HvDRR33,ind=c(1:nrow(HvDRR33$pheno))[!is.na(HvDRR33$pheno)])

HvDRR33per <- calc.genoprob(HvDRR33, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR33out <- scanone(HvDRR33per, method = "hk")
plot(HvDRR33out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR33permout <- scanone(HvDRR33per, method = "hk", n.perm = 4000)
summary(HvDRR33out, perms = HvDRR33permout, alpha = 0.05, pvalues = TRUE)
HvDRR33qtl1 <- makeqtl(HvDRR33per, chr = 2, pos = 95, what = "prob")
HvDRR33outc2 <- addqtl(HvDRR33per, qtl = HvDRR33qtl1, method = "hk")
summary(HvDRR33outc2, perms = HvDRR33permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR33out, HvDRR33outc2, chr = -2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR33out, HvDRR33outc2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)

HvDRR33rqtl <- refineqtl(HvDRR33per, qtl = HvDRR33qtl1, method = "hk", verbose = FALSE)
summary(HvDRR33rqtl)
summary(fitqtl(HvDRR33per, qtl = HvDRR33rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR33per, 2, 95)
effectplot(HvDRR33per, mname1 = "JHI-Hv50k-2016-75950")
summary(HvDRR33per, mname1 = "JHI-Hv50k-2016-75950")
plotPXG(HvDRR33per, "JHI-Hv50k-2016-75950")

summary(fitqtl(HvDRR33per, qtl = HvDRR33rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR33CI2 <- lodint(HvDRR33rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(-6.15, 78, lower.tail = FALSE)
#p<0.01


#HvDRR34####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR34 <- Populations$HvDRR34

HvDRR34$geno$'1'$data<- apply(HvDRR34$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'2'$data<- apply(HvDRR34$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'3'$data<- apply(HvDRR34$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'4'$data<- apply(HvDRR34$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'5'$data<- apply(HvDRR34$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'6'$data<- apply(HvDRR34$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'7'$data<- apply(HvDRR34$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR34$geno$'1'$data<- apply(HvDRR34$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'2'$data<- apply(HvDRR34$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'3'$data<- apply(HvDRR34$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'4'$data<- apply(HvDRR34$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'5'$data<- apply(HvDRR34$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'6'$data<- apply(HvDRR34$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR34$geno$'7'$data<- apply(HvDRR34$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR34)[1] <- "riself"

summary(HvDRR34)
plot(HvDRR34)
HvDRR34 <- subset(HvDRR34,ind=c(1:nrow(HvDRR34$pheno))[!is.na(HvDRR34$pheno)])

HvDRR34per <- calc.genoprob(HvDRR34, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR34out <- scanone(HvDRR34per, method = "hk")
plot(HvDRR34out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR34permout <- scanone(HvDRR34per, method = "hk", n.perm=4000)
summary(HvDRR34out, perms = HvDRR34permout, alpha = 0.05, pvalues = TRUE)
HvDRR34qtl1 <- makeqtl(HvDRR34per, chr = 2, pos = 250, what = "prob")
HvDRR34outc2 <- addqtl(HvDRR34per, qtl = HvDRR34qtl1, method = "hk")
summary(HvDRR34outc2, perms = HvDRR34permout, alpha = 0.05, pvalues = TRUE)

HvDRR34rqtl <- refineqtl(HvDRR34per, qtl = HvDRR34qtl1, method = "hk", verbose = FALSE)
summary(HvDRR34rqtl)
summary(fitqtl(HvDRR34per, qtl = HvDRR34rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR34per, 2, 250.2)
effectplot(HvDRR34per, mname1 = "SCRI_RS_55841")
summary(HvDRR34per, mname1 = "SCRI_RS_55841")
plotPXG(HvDRR34per, "SCRI_RS_55841")

summary(fitqtl(HvDRR34per, qtl = HvDRR34rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR34CI2 <- lodint(HvDRR34rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(4.827, 32, lower.tail = FALSE)
#p<0.01

#HvDRR35####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR35 <- Populations$HvDRR35

HvDRR35$geno$'1'$data<- apply(HvDRR35$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'2'$data<- apply(HvDRR35$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'3'$data<- apply(HvDRR35$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'4'$data<- apply(HvDRR35$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'5'$data<- apply(HvDRR35$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'6'$data<- apply(HvDRR35$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'7'$data<- apply(HvDRR35$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR35$geno$'1'$data<- apply(HvDRR35$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'2'$data<- apply(HvDRR35$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'3'$data<- apply(HvDRR35$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'4'$data<- apply(HvDRR35$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'5'$data<- apply(HvDRR35$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'6'$data<- apply(HvDRR35$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR35$geno$'7'$data<- apply(HvDRR35$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR35)[1] <- "riself"


summary(HvDRR35)
plot(HvDRR35)
HvDRR35 <- subset(HvDRR35,ind=c(1:nrow(HvDRR35$pheno))[!is.na(HvDRR35$pheno)])
HvDRR35per <- calc.genoprob(HvDRR35, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR35out <- scanone(HvDRR35per, method = "hk")
plot(HvDRR35out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR35permout <- scanone(HvDRR35per, method = "hk", n.perm=4000)
summary(HvDRR35out, perms = HvDRR35permout, alpha = 0.05, pvalues = TRUE)

#HvDRR36####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR36 <- Populations$HvDRR36

HvDRR36$geno$'1'$data<- apply(HvDRR36$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'2'$data<- apply(HvDRR36$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'3'$data<- apply(HvDRR36$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'4'$data<- apply(HvDRR36$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'5'$data<- apply(HvDRR36$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'6'$data<- apply(HvDRR36$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'7'$data<- apply(HvDRR36$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR36$geno$'1'$data<- apply(HvDRR36$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'2'$data<- apply(HvDRR36$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'3'$data<- apply(HvDRR36$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'4'$data<- apply(HvDRR36$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'5'$data<- apply(HvDRR36$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'6'$data<- apply(HvDRR36$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR36$geno$'7'$data<- apply(HvDRR36$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR36)[1] <- "riself"

summary(HvDRR36)
plot(HvDRR36)
HvDRR36 <- subset(HvDRR36,ind=c(1:nrow(HvDRR36$pheno))[!is.na(HvDRR36$pheno)])

HvDRR36per <- calc.genoprob(HvDRR36, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR36out <- scanone(HvDRR36per, method = "hk")
plot(HvDRR36out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR36permout <- scanone(HvDRR36per, method = "hk", n.perm=4000)
summary(HvDRR36out, perms = HvDRR36permout, alpha = 0.05, pvalues = TRUE)

#HvDRR37####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR37 <- Populations$HvDRR37

HvDRR37$geno$'1'$data<- apply(HvDRR37$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'2'$data<- apply(HvDRR37$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'3'$data<- apply(HvDRR37$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'4'$data<- apply(HvDRR37$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'5'$data<- apply(HvDRR37$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'6'$data<- apply(HvDRR37$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'7'$data<- apply(HvDRR37$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR37$geno$'1'$data<- apply(HvDRR37$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'2'$data<- apply(HvDRR37$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'3'$data<- apply(HvDRR37$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'4'$data<- apply(HvDRR37$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'5'$data<- apply(HvDRR37$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'6'$data<- apply(HvDRR37$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR37$geno$'7'$data<- apply(HvDRR37$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR37)[1] <- "riself"

summary(HvDRR37)
plot(HvDRR37)
HvDRR37 <- subset(HvDRR37,ind=c(1:nrow(HvDRR37$pheno))[!is.na(HvDRR37$pheno)])

HvDRR37per <- calc.genoprob(HvDRR37, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR37out <- scanone(HvDRR37per, method = "hk")
plot(HvDRR37out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR37permout <- scanone(HvDRR37per, method = "hk", n.perm=4000)
summary(HvDRR37out, perms = HvDRR37permout, alpha = 0.05, pvalues = TRUE)
HvDRR37qtl1 <- makeqtl(HvDRR37per, chr = 2, pos = 286, what = "prob")
HvDRR37outc2 <- addqtl(HvDRR37per, qtl = HvDRR37qtl1, method = "hk")
summary(HvDRR37outc2, perms = HvDRR37permout, alpha = 0.05, pvalues = TRUE)
HvDRR37qtl2 <- makeqtl(HvDRR37per, chr = c(2,6), pos = c(286,110), what = "prob")
HvDRR37outc26 <- addqtl(HvDRR37per, qtl = HvDRR37qtl2, method = "hk")
summary(HvDRR37outc26, perms = HvDRR37permout, alpha = 0.05, pvalues = TRUE)

HvDRR37rqtl <- refineqtl(HvDRR37per, qtl = HvDRR37qtl2, method = "hk", verbose = FALSE)
summary(HvDRR37rqtl)
summary(fitqtl(HvDRR37per, qtl = HvDRR37rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR37per, 2, 286)
effectplot(HvDRR37per, mname1 = "JHI-Hv50k-2016-131048")
summary(HvDRR37per, mname1 = "JHI-Hv50k-2016-131048")
plotPXG(HvDRR37per, "JHI-Hv50k-2016-131048")

find.marker(HvDRR37per, 6, 110)
effectplot(HvDRR37per, mname1 = "JHI-Hv50k-2016-401443")
summary(HvDRR37per, mname1 = "JHI-Hv50k-2016-401443")
plotPXG(HvDRR37per, "JHI-Hv50k-2016-401443")

summary(fitqtl(HvDRR37per, qtl = HvDRR37rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2))

HvDRR37CI2 <- lodint(HvDRR37rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR37CI6 <- lodint(HvDRR37rqtl, chr = 6, qtl.index = 2, expandtomarkers = TRUE)

pt(5.654, 60, lower.tail = FALSE)
pt(-4.462, 60, lower.tail = FALSE)

#HvDRR38####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR38 <- Populations$HvDRR38

HvDRR38$geno$'1'$data<- apply(HvDRR38$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'2'$data<- apply(HvDRR38$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'3'$data<- apply(HvDRR38$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'4'$data<- apply(HvDRR38$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'5'$data<- apply(HvDRR38$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'6'$data<- apply(HvDRR38$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'7'$data<- apply(HvDRR38$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR38$geno$'1'$data<- apply(HvDRR38$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'2'$data<- apply(HvDRR38$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'3'$data<- apply(HvDRR38$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'4'$data<- apply(HvDRR38$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'5'$data<- apply(HvDRR38$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'6'$data<- apply(HvDRR38$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR38$geno$'7'$data<- apply(HvDRR38$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR38)[1] <- "riself"

summary(HvDRR38)
plot(HvDRR38)
HvDRR38 <- subset(HvDRR38,ind=c(1:nrow(HvDRR38$pheno))[!is.na(HvDRR38$pheno)])

HvDRR38per <- calc.genoprob(HvDRR38, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR38out <- scanone(HvDRR38per, method = "hk")
plot(HvDRR38out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR38permout <- scanone(HvDRR38per, method = "hk", n.perm=4000)
summary(HvDRR38out, perms = HvDRR38permout, alpha = 0.05, pvalues = TRUE)
HvDRR38qtl1 <- makeqtl(HvDRR38per, chr = 2, pos = 363, what = "prob")
HvDRR38outc2 <- addqtl(HvDRR38per, qtl = HvDRR38qtl1, method = "hk")
summary(HvDRR38outc2, perms = HvDRR38permout, alpha = 0.05, pvalues = TRUE)

HvDRR38rqtl <- refineqtl(HvDRR38per, qtl = HvDRR38qtl1, method = "hk", verbose = FALSE)
summary(HvDRR38rqtl)
summary(fitqtl(HvDRR38per, qtl = HvDRR38rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR38per, 2, 363)
effectplot(HvDRR38per, mname1 = "SCRI_RS_147230")
summary(HvDRR38per, mname1 = "SCRI_RS_147230")
plotPXG(HvDRR38per, "SCRI_RS_147230")

summary(fitqtl(HvDRR38per, qtl = HvDRR38rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR38CI2 <- lodint(HvDRR38rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)

pt(7.845, 74, lower.tail = FALSE)
#p<0.01

#HvDRR39####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR39 <- Populations$HvDRR39

HvDRR39$geno$'1'$data<- apply(HvDRR39$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'2'$data<- apply(HvDRR39$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'3'$data<- apply(HvDRR39$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'4'$data<- apply(HvDRR39$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'5'$data<- apply(HvDRR39$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'6'$data<- apply(HvDRR39$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'7'$data<- apply(HvDRR39$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR39$geno$'1'$data<- apply(HvDRR39$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'2'$data<- apply(HvDRR39$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'3'$data<- apply(HvDRR39$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'4'$data<- apply(HvDRR39$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'5'$data<- apply(HvDRR39$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'6'$data<- apply(HvDRR39$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR39$geno$'7'$data<- apply(HvDRR39$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR39)[1] <- "riself"

summary(HvDRR39)
plot(HvDRR39)
HvDRR39 <- subset(HvDRR39,ind=c(1:nrow(HvDRR39$pheno))[!is.na(HvDRR39$pheno)])

testmarkername <- as.data.frame(HvDRR39$geno$`2`$map)
which(row.names(testmarkername) == 'JHI-Hv50k-2016-73673')
marker2drop <- markernames(HvDRR39, chr = 2)[54]
HvDRR39 <- drop.markers(HvDRR39, marker2drop)

HvDRR39per <- calc.genoprob(HvDRR39, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR39out <- scanone(HvDRR39per, method = "hk")
plot(HvDRR39out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR39permout <- scanone(HvDRR39per, method = "hk", n.perm=4000)
summary(HvDRR39out, perms = HvDRR39permout, alpha = 0.05, pvalues = TRUE)
HvDRR39qtl1 <- makeqtl(HvDRR39per, chr = 2, pos = 239, what = "prob")
HvDRR39outc2 <- addqtl(HvDRR39per, qtl = HvDRR39qtl1, method = "hk")
summary(HvDRR39outc2, perms = HvDRR39permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR39out, HvDRR39outc2, chr = -2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
plot(HvDRR39out, HvDRR39outc2, col = c("blue", "red"), ylab = "LOD score", ylim = c(0,5), alternate.chrid = TRUE)
HvDRR39qtl2 <- makeqtl(HvDRR39per, c(2,2), c(239,41.6), what = "prob")
HvDRR39outc22 <- addqtl(HvDRR39per, qtl = HvDRR39qtl2, method = "hk")
summary(HvDRR39outc22, perms = HvDRR39permout, alpha = 0.05, pvalues = TRUE)

HvDRR39rqtl <- refineqtl(HvDRR39per, qtl = HvDRR39qtl2, method = "hk", verbose = FALSE)
summary(HvDRR39rqtl)
summary(fitqtl(HvDRR39per, qtl = HvDRR39rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR39per, 2, 238.8)
effectplot(HvDRR39per, mname1 = "JHI-Hv50k-2016-130184")
summary(HvDRR39per, mname1 = "JHI-Hv50k-2016-130184")
plotPXG(HvDRR39per, "JHI-Hv50k-2016-130184")

find.marker(HvDRR39per, 2, 41.6)
effectplot(HvDRR39per, mname1 = "JHI-Hv50k-2016-73618")
summary(HvDRR39per, mname1 = "JHI-Hv50k-2016-73618")
plotPXG(HvDRR39per, "JHI-Hv50k-2016-73618")

summary(fitqtl(HvDRR39per, qtl = HvDRR39rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2))

HvDRR39CI2 <- lodint(HvDRR39rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR39CI22 <- lodint(HvDRR39rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)

pt(9.362, 102, lower.tail = FALSE)
#p<0.01
pt(4.726, 102, lower.tail = FALSE)
#p<0.01

#HvDRR40####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR40 <- Populations$HvDRR40

HvDRR40$geno$'1'$data<- apply(HvDRR40$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'2'$data<- apply(HvDRR40$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'3'$data<- apply(HvDRR40$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'4'$data<- apply(HvDRR40$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'5'$data<- apply(HvDRR40$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'6'$data<- apply(HvDRR40$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'7'$data<- apply(HvDRR40$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR40$geno$'1'$data<- apply(HvDRR40$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'2'$data<- apply(HvDRR40$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'3'$data<- apply(HvDRR40$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'4'$data<- apply(HvDRR40$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'5'$data<- apply(HvDRR40$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'6'$data<- apply(HvDRR40$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR40$geno$'7'$data<- apply(HvDRR40$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR40)[1] <- "riself"

summary(HvDRR40)
plot(HvDRR40)
HvDRR40 <- subset(HvDRR40,ind=c(1:nrow(HvDRR40$pheno))[!is.na(HvDRR40$pheno)])

HvDRR40per <- calc.genoprob(HvDRR40, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR40out <- scanone(HvDRR40per, method = "hk")
plot(HvDRR40out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR40permout <- scanone(HvDRR40per, method = "hk", n.perm=4000)
summary(HvDRR40out, perms = HvDRR40permout, alpha = 0.05, pvalues = TRUE)

#HvDRR41####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR41 <- Populations$HvDRR41

HvDRR41$geno$'1'$data<- apply(HvDRR41$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'2'$data<- apply(HvDRR41$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'3'$data<- apply(HvDRR41$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'4'$data<- apply(HvDRR41$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'5'$data<- apply(HvDRR41$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'6'$data<- apply(HvDRR41$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'7'$data<- apply(HvDRR41$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR41$geno$'1'$data<- apply(HvDRR41$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'2'$data<- apply(HvDRR41$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'3'$data<- apply(HvDRR41$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'4'$data<- apply(HvDRR41$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'5'$data<- apply(HvDRR41$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'6'$data<- apply(HvDRR41$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR41$geno$'7'$data<- apply(HvDRR41$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR41)[1] <- "riself"

summary(HvDRR41)
plot(HvDRR41)
HvDRR41 <- subset(HvDRR41,ind=c(1:nrow(HvDRR41$pheno))[!is.na(HvDRR41$pheno)])

HvDRR41per <- calc.genoprob(HvDRR41, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR41out <- scanone(HvDRR41per, method = "hk")
plot(HvDRR41out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR41permout <- scanone(HvDRR41per, method = "hk", n.perm=4000)
summary(HvDRR41out, perms = HvDRR41permout, alpha = 0.05, pvalues = TRUE)
HvDRR41qtl1 <- makeqtl(HvDRR41per, chr = c(1,2), pos = c(136,45.8), what = "prob")
HvDRR41outc12 <- addqtl(HvDRR41per, qtl = HvDRR41qtl1, method = "hk")
summary(HvDRR41outc12, perms = HvDRR41permout, alpha = 0.05, pvalues = TRUE)
HvDRR41qtl2 <- makeqtl(HvDRR41per, chr = c(1,2,2,7), pos = c(136,45.8,151.9,47.8), what = "prob")
HvDRR41outc1227 <- addqtl(HvDRR41per, qtl = HvDRR41qtl2, method = "hk")
summary(HvDRR41outc1227, perms = HvDRR41permout, alpha = 0.05, pvalues = TRUE)
HvDRR41qtl3 <- makeqtl(HvDRR41per, chr = c(1,2,2,7,2,7), pos = c(136,45.8,151.9,47.8,116,97.9), what = "prob")
HvDRR41outc122727 <- addqtl(HvDRR41per, qtl = HvDRR41qtl3, method = "hk")
summary(HvDRR41outc122727, perms = HvDRR41permout, alpha = 0.05, pvalues = TRUE)

HvDRR41rqtl <- refineqtl(HvDRR41per, qtl = HvDRR41qtl3, method = "hk", verbose = FALSE)
summary(HvDRR41rqtl)
summary(fitqtl(HvDRR41per, qtl = HvDRR41rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR41per, 1, 147.4)
effectplot(HvDRR41per, mname1 = "BOPA2_12_30191")
summary(HvDRR41per, mname1 = "BOPA2_12_30191")
plotPXG(HvDRR41per, "BOPA2_12_30191")

find.marker(HvDRR41per, 2, 45.8)
effectplot(HvDRR41per, mname1 = "JHI-Hv50k-2016-73417")
summary(HvDRR41per, mname1 = "JHI-Hv50k-2016-73417")
plotPXG(HvDRR41per, "JHI-Hv50k-2016-73417")

find.marker(HvDRR41per, 2, 148.6)
effectplot(HvDRR41per, mname1 = "JHI-Hv50k-2016-102289")
summary(HvDRR41per, mname1 = "JHI-Hv50k-2016-102289")
plotPXG(HvDRR41per, "JHI-Hv50k-2016-102289")

find.marker(HvDRR41per, 7, 51)
effectplot(HvDRR41per, mname1 = "JHI-Hv50k-2016-460797")
summary(HvDRR41per, mname1 = "JHI-Hv50k-2016-460797")
plotPXG(HvDRR41per, "JHI-Hv50k-2016-460797")

find.marker(HvDRR41per, 2, 115.2)
effectplot(HvDRR41per, mname1 = "BOPA2_12_10330")
summary(HvDRR41per, mname1 = "BOPA2_12_10330")
plotPXG(HvDRR41per, "BOPA2_12_10330")

find.marker(HvDRR41per, 7, 97.9)
effectplot(HvDRR41per, mname1 = "JHI-Hv50k-2016-471764")
summary(HvDRR41per, mname1 = "JHI-Hv50k-2016-471764")
plotPXG(HvDRR41per, "JHI-Hv50k-2016-471764")

summary(fitqtl(HvDRR41per, qtl = HvDRR41rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6))

HvDRR41CI1 <- lodint(HvDRR41rqtl, chr = 1, qtl.index = 1, expandtomarkers = TRUE)
HvDRR41CI2 <- lodint(HvDRR41rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)
HvDRR41CI22 <- lodint(HvDRR41rqtl, chr = 2, qtl.index = 3, expandtomarkers = TRUE)
HvDRR41CI7 <- lodint(HvDRR41rqtl, chr = 7, qtl.index = 4, expandtomarkers = TRUE)
HvDRR41CI222 <- lodint(HvDRR41rqtl, chr = 2, qtl.index = 5, expandtomarkers = TRUE)
HvDRR41CI77 <- lodint(HvDRR41rqtl, chr = 7, qtl.index = 6, expandtomarkers = TRUE)

pt(3.994, 86, lower.tail = FALSE)
pt(11.343, 86, lower.tail = FALSE)
pt(-5.401, 86, lower.tail = FALSE)
pt(7.825, 86, lower.tail = FALSE)
pt(4.210, 86, lower.tail = FALSE)
pt(-5.469, 86, lower.tail = FALSE)

#HvDRR42####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR42 <- Populations$HvDRR42

HvDRR42$geno$'1'$data<- apply(HvDRR42$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'2'$data<- apply(HvDRR42$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'3'$data<- apply(HvDRR42$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'4'$data<- apply(HvDRR42$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'5'$data<- apply(HvDRR42$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'6'$data<- apply(HvDRR42$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'7'$data<- apply(HvDRR42$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR42$geno$'1'$data<- apply(HvDRR42$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'2'$data<- apply(HvDRR42$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'3'$data<- apply(HvDRR42$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'4'$data<- apply(HvDRR42$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'5'$data<- apply(HvDRR42$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'6'$data<- apply(HvDRR42$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR42$geno$'7'$data<- apply(HvDRR42$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR42)[1] <- "riself"

summary(HvDRR42)
plot(HvDRR42)
HvDRR42 <- subset(HvDRR42,ind=c(1:nrow(HvDRR42$pheno))[!is.na(HvDRR42$pheno)])

HvDRR42per <- calc.genoprob(HvDRR42, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR42out <- scanone(HvDRR42per, method = "hk")
plot(HvDRR42out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR42permout <- scanone(HvDRR42per, method = "hk", n.perm=4000)
summary(HvDRR42out, perms = HvDRR42permout, alpha = 0.05, pvalues = TRUE)

#HvDRR43####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR43 <- Populations$HvDRR43

HvDRR43$geno$'1'$data<- apply(HvDRR43$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'2'$data<- apply(HvDRR43$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'3'$data<- apply(HvDRR43$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'4'$data<- apply(HvDRR43$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'5'$data<- apply(HvDRR43$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'6'$data<- apply(HvDRR43$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'7'$data<- apply(HvDRR43$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR43$geno$'1'$data<- apply(HvDRR43$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'2'$data<- apply(HvDRR43$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'3'$data<- apply(HvDRR43$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'4'$data<- apply(HvDRR43$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'5'$data<- apply(HvDRR43$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'6'$data<- apply(HvDRR43$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR43$geno$'7'$data<- apply(HvDRR43$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR43)[1] <- "riself"

summary(HvDRR43)
plot(HvDRR43)
HvDRR43 <- subset(HvDRR43,ind=c(1:nrow(HvDRR43$pheno))[!is.na(HvDRR43$pheno)])

HvDRR43per <- calc.genoprob(HvDRR43, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR43out <- scanone(HvDRR43per, method = "hk")
plot(HvDRR43out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR43permout <- scanone(HvDRR43per, method = "hk", n.perm=4000)
summary(HvDRR43out, perms = HvDRR43permout, alpha = 0.05, pvalues = TRUE)
HvDRR43qtl1 <- makeqtl(HvDRR43per, chr = 4, pos = 131, what = "prob")
HvDRR43outc4 <- addqtl(HvDRR43per, qtl = HvDRR43qtl1, method = "hk")
summary(HvDRR43outc4, perms = HvDRR43permout, alpha = 0.05, pvalues = TRUE)

HvDRR43rqtl <- refineqtl(HvDRR43per, qtl = HvDRR43qtl1, method = "hk", verbose = FALSE)
summary(HvDRR43rqtl)
summary(fitqtl(HvDRR43per, qtl = HvDRR43rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR43per, 4, 131.1)
effectplot(HvDRR43per, mname1 = "JHI-Hv50k-2016-268870")
summary(HvDRR43per, mname1 = "JHI-Hv50k-2016-268870")
plotPXG(HvDRR43per, "JHI-Hv50k-2016-268870")

summary(fitqtl(HvDRR43per, qtl = HvDRR43rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR43CI <- lodint(HvDRR43rqtl, chr = 4, qtl.index = 1, expandtomarkers = TRUE)

pt(5.239, 125, lower.tail = FALSE)

#HvDRR44####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR44 <- Populations$HvDRR44

HvDRR44$geno$'1'$data<- apply(HvDRR44$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'2'$data<- apply(HvDRR44$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'3'$data<- apply(HvDRR44$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'4'$data<- apply(HvDRR44$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'5'$data<- apply(HvDRR44$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'6'$data<- apply(HvDRR44$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'7'$data<- apply(HvDRR44$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR44$geno$'1'$data<- apply(HvDRR44$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'2'$data<- apply(HvDRR44$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'3'$data<- apply(HvDRR44$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'4'$data<- apply(HvDRR44$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'5'$data<- apply(HvDRR44$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'6'$data<- apply(HvDRR44$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR44$geno$'7'$data<- apply(HvDRR44$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR44)[1] <- "riself"

summary(HvDRR44)
plot(HvDRR44)
HvDRR44 <- subset(HvDRR44,ind=c(1:nrow(HvDRR44$pheno))[!is.na(HvDRR44$pheno)])

HvDRR44per <- calc.genoprob(HvDRR44, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR44out <- scanone(HvDRR44per, method = "hk")
plot(HvDRR44out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR44permout <- scanone(HvDRR44per, method = "hk", n.perm = 4000)
summary(HvDRR44out, perms = HvDRR44permout, alpha = 0.05, pvalues = TRUE)
HvDRR44qtl1 <- makeqtl(HvDRR44per, chr = 5, pos = 262, what = "prob")
HvDRR44outc5 <- addqtl(HvDRR44per, qtl = HvDRR44qtl1, method = "hk")
summary(HvDRR44outc5, perms = HvDRR44permout, alpha = 0.05, pvalues = TRUE)

HvDRR44rqtl <- refineqtl(HvDRR44per, qtl = HvDRR44qtl1, method = "hk", verbose = FALSE)
summary(HvDRR44rqtl)
summary(fitqtl(HvDRR44per, qtl = HvDRR44rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR44per, 5, 262)
effectplot(HvDRR44per, mname1 = "JHI-Hv50k-2016-338130")
summary(HvDRR44per, mname1 = "JHI-Hv50k-2016-338130")
plotPXG(HvDRR44per, "JHI-Hv50k-2016-338130")

summary(fitqtl(HvDRR44per, qtl = HvDRR44rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1))

HvDRR44CI <- lodint(HvDRR44rqtl, chr = 5, qtl.index = 1, expandtomarkers = TRUE)

pt(-6.108, 95, lower.tail = FALSE)
#p<0.01

#HvDRR45####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR45 <- Populations$HvDRR45

HvDRR45$geno$'1'$data<- apply(HvDRR45$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'2'$data<- apply(HvDRR45$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'3'$data<- apply(HvDRR45$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'4'$data<- apply(HvDRR45$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'5'$data<- apply(HvDRR45$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'6'$data<- apply(HvDRR45$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'7'$data<- apply(HvDRR45$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR45$geno$'1'$data<- apply(HvDRR45$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'2'$data<- apply(HvDRR45$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'3'$data<- apply(HvDRR45$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'4'$data<- apply(HvDRR45$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'5'$data<- apply(HvDRR45$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'6'$data<- apply(HvDRR45$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR45$geno$'7'$data<- apply(HvDRR45$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR45)[1] <- "riself"

summary(HvDRR45)
plot(HvDRR45)
HvDRR45 <- subset(HvDRR45,ind=c(1:nrow(HvDRR45$pheno))[!is.na(HvDRR45$pheno)])

HvDRR45per <- calc.genoprob(HvDRR45, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR45out <- scanone(HvDRR45per, method = "hk")
plot(HvDRR45out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR45permout <- scanone(HvDRR45per, method = "hk", n.perm=4000)
summary(HvDRR45out, perms = HvDRR45permout, alpha = 0.05, pvalues = TRUE)
HvDRR45qtl1 <- makeqtl(HvDRR45per, chr = 5, pos = 196, what = "prob")
HvDRR45outc5 <- addqtl(HvDRR45per, qtl = HvDRR45qtl1, method = "hk")
summary(HvDRR45outc5, perms = HvDRR45permout, alpha = 0.05, pvalues = TRUE)
HvDRR45qtl2 <- makeqtl(HvDRR45per, chr = c(5,1), pos = c(196,90.7), what = "prob")
HvDRR45outc51 <- addqtl(HvDRR45per, qtl = HvDRR45qtl2, method = "hk")
summary(HvDRR45outc51, perms = HvDRR45permout, alpha = 0.05, pvalues = TRUE)
HvDRR45qtl3 <- makeqtl(HvDRR45per, chr = c(5,1,7), pos = c(196,90.7,140), what = "prob")
HvDRR45outc517 <- addqtl(HvDRR45per, qtl = HvDRR45qtl3, method = "hk")
summary(HvDRR45outc517, perms = HvDRR45permout, alpha = 0.05, pvalues = TRUE)

HvDRR45rqtl <- refineqtl(HvDRR45per, qtl = HvDRR45qtl3, method = "hk", verbose = FALSE)
summary(HvDRR45rqtl)
summary(fitqtl(HvDRR45per, qtl = HvDRR45rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR45per, 5, 196)
effectplot(HvDRR45per, mname1 = "JHI-Hv50k-2016-333733")
summary(HvDRR45per, mname1 = "JHI-Hv50k-2016-333733")
plotPXG(HvDRR45per, "JHI-Hv50k-2016-333733")

find.marker(HvDRR45per, 1, 90.7)
effectplot(HvDRR45per, mname1 = "JHI-Hv50k-2016-20787")
summary(HvDRR45per, mname1 = "JHI-Hv50k-2016-20787")
plotPXG(HvDRR45per, "JHI-Hv50k-2016-20787")

find.marker(HvDRR45per, 7, 139.5)
effectplot(HvDRR45per, mname1 = "JHI-Hv50k-2016-491516")
summary(HvDRR45per, mname1 = "JHI-Hv50k-2016-491516")
plotPXG(HvDRR45per, "JHI-Hv50k-2016-491516")

summary(fitqtl(HvDRR45per, qtl = HvDRR45rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2 + Q3))

HvDRR45CI5 <- lodint(HvDRR45rqtl, chr = 5, qtl.index = 1, expandtomarkers = TRUE)
HvDRR45CI1 <- lodint(HvDRR45rqtl, chr = 1, qtl.index = 2, expandtomarkers = TRUE)
HvDRR45CI7 <- lodint(HvDRR45rqtl, chr = 7, qtl.index = 3, expandtomarkers = TRUE)

pt(-7.547, 73, lower.tail = FALSE)
#p<0.01
pt(4.586, 73, lower.tail = FALSE)
#p<0.01
pt(4.194, 73, lower.tail = FALSE)
#p<0.01

#HvDRR46####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR46 <- Populations$HvDRR46

HvDRR46$geno$'1'$data<- apply(HvDRR46$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'2'$data<- apply(HvDRR46$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'3'$data<- apply(HvDRR46$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'4'$data<- apply(HvDRR46$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'5'$data<- apply(HvDRR46$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'6'$data<- apply(HvDRR46$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'7'$data<- apply(HvDRR46$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR46$geno$'1'$data<- apply(HvDRR46$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'2'$data<- apply(HvDRR46$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'3'$data<- apply(HvDRR46$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'4'$data<- apply(HvDRR46$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'5'$data<- apply(HvDRR46$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'6'$data<- apply(HvDRR46$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR46$geno$'7'$data<- apply(HvDRR46$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR46)[1] <- "riself"

summary(HvDRR46)
plot(HvDRR46)
HvDRR46 <- subset(HvDRR46,ind=c(1:nrow(HvDRR46$pheno))[!is.na(HvDRR46$pheno)])

HvDRR46per <- calc.genoprob(HvDRR46, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR46out <- scanone(HvDRR46per, method = "hk")
plot(HvDRR46out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR46permout <- scanone(HvDRR46per, method = "hk", n.perm=4000)
summary(HvDRR46out, perms = HvDRR46permout, alpha = 0.05, pvalues = TRUE)

#HvDRR47####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR47 <- Populations$HvDRR47

HvDRR47$geno$'1'$data<- apply(HvDRR47$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'2'$data<- apply(HvDRR47$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'3'$data<- apply(HvDRR47$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'4'$data<- apply(HvDRR47$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'5'$data<- apply(HvDRR47$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'6'$data<- apply(HvDRR47$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'7'$data<- apply(HvDRR47$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR47$geno$'1'$data<- apply(HvDRR47$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'2'$data<- apply(HvDRR47$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'3'$data<- apply(HvDRR47$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'4'$data<- apply(HvDRR47$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'5'$data<- apply(HvDRR47$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'6'$data<- apply(HvDRR47$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR47$geno$'7'$data<- apply(HvDRR47$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR47)[1] <- "riself"

summary(HvDRR47)
plot(HvDRR47)
HvDRR47 <- subset(HvDRR47,ind=c(1:nrow(HvDRR47$pheno))[!is.na(HvDRR47$pheno)])

testmarkername <- as.data.frame(HvDRR47$geno$`3`$map)
which(row.names(testmarkername) == 'SCRI_RS_106728')
marker2drop <- markernames(HvDRR47, chr = 3)[532]
HvDRR47 <- drop.markers(HvDRR47, marker2drop)
HvDRR47per <- calc.genoprob(HvDRR47, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR47out <- scanone(HvDRR47per, method = "hk")
plot(HvDRR47out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR47permout <- scanone(HvDRR47per, method = "hk", n.perm = 4000)
summary(HvDRR47out, perms = HvDRR47permout, alpha = 0.05, pvalues = TRUE)
HvDRR47qtl1 <- makeqtl(HvDRR47per, 2, 76, what = "prob")
HvDRR47outc2 <- addqtl(HvDRR47per, qtl = HvDRR47qtl1, method = "hk")
summary(HvDRR47outc2, perms = HvDRR47permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR47out, HvDRR47outc2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)
HvDRR47qtl2 <- makeqtl(HvDRR47per, c(2,3), c(76,340), what = "prob")
HvDRR47outc23 <- addqtl(HvDRR47per, qtl = HvDRR47qtl2, method = "hk")
summary(HvDRR47outc23, perms = HvDRR47permout, alpha = 0.05, pvalues = TRUE)
plot(HvDRR47out, HvDRR47outc2, col = c("blue", "red"), ylab = "LOD score", alternate.chrid = TRUE)

HvDRR47rqtl <- refineqtl(HvDRR47per, qtl = HvDRR47qtl2, method = "hk", verbose = FALSE)
summary(HvDRR47rqtl)
summary(fitqtl(HvDRR47per, qtl = HvDRR47rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR47per, 2, 74.8)
effectplot(HvDRR47per, mname1 = "JHI-Hv50k-2016-77567")
summary(HvDRR47per, mname1 = "JHI-Hv50k-2016-77567")
plotPXG(HvDRR47per, "JHI-Hv50k-2016-77567")

find.marker(HvDRR47per, 3, 339.8)
effectplot(HvDRR47per, mname1 = "JHI-Hv50k-2016-211559")
summary(HvDRR47per, mname1 = "JHI-Hv50k-2016-211559")
plotPXG(HvDRR47per, "JHI-Hv50k-2016-211559")

summary(fitqtl(HvDRR47per, qtl = HvDRR47rqtl, method = "hk", get.ests=TRUE, formula=Flowering~ Q1 + Q2))

HvDRR47CI2 <- lodint(HvDRR47rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR47CI3 <- lodint(HvDRR47rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)

pt(5.804, 97, lower.tail = FALSE)
#p<0.01
pt(-2.8329, 97, lower.tail = FALSE)
#not significant

#HvDRR48####

Populations <- readRDS("Population_analysis_Plant_Height.RDS")
library(qtl)
HvDRR48 <- Populations$HvDRR48

HvDRR48$geno$'1'$data<- apply(HvDRR48$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'2'$data<- apply(HvDRR48$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'3'$data<- apply(HvDRR48$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'4'$data<- apply(HvDRR48$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'5'$data<- apply(HvDRR48$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'6'$data<- apply(HvDRR48$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'7'$data<- apply(HvDRR48$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
HvDRR48$geno$'1'$data<- apply(HvDRR48$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'2'$data<- apply(HvDRR48$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'3'$data<- apply(HvDRR48$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'4'$data<- apply(HvDRR48$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'5'$data<- apply(HvDRR48$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'6'$data<- apply(HvDRR48$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
HvDRR48$geno$'7'$data<- apply(HvDRR48$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})

class(HvDRR48)[1] <- "riself"


summary(HvDRR48)
plot(HvDRR48)
HvDRR48 <- subset(HvDRR48,ind=c(1:nrow(HvDRR48$pheno))[!is.na(HvDRR48$pheno)])

HvDRR48per <- calc.genoprob(HvDRR48, step = 1, error.prob = 0.01, map.function = "haldane")
HvDRR48out <- scanone(HvDRR48per, method = "hk")
plot(HvDRR48out, ylab = "LOD score", alternate.chrid = TRUE)
HvDRR48permout <- scanone(HvDRR48per, method = "hk", n.perm=4000)
summary(HvDRR48out, perms = HvDRR48permout, alpha = 0.05, pvalues = TRUE)
HvDRR48qtl1 <- makeqtl(HvDRR48per, chr = 2, pos = 40.5, what = "prob")
HvDRR48outc2 <- addqtl(HvDRR48per, qtl = HvDRR48qtl1, method = "hk")
summary(HvDRR48outc2, perms = HvDRR48permout, alpha = 0.05, pvalues = TRUE)
HvDRR48qtl2 <- makeqtl(HvDRR48per, chr = c(2,4), pos = c(40.5,41.7), what = "prob")
HvDRR48outc24 <- addqtl(HvDRR48per, qtl = HvDRR48qtl2, method = "hk")
summary(HvDRR48outc24, perms = HvDRR48permout, alpha = 0.05, pvalues = TRUE)

HvDRR48rqtl <- refineqtl(HvDRR48per, qtl = HvDRR48qtl2, method = "hk", verbose = FALSE)
summary(HvDRR48rqtl)
summary(fitqtl(HvDRR48per, qtl = HvDRR48rqtl, method = "hk"), pvalues = FALSE)

find.marker(HvDRR48per, 2, 35.9)
effectplot(HvDRR48per, mname1 = "JHI-Hv50k-2016-73692")
summary(HvDRR48per, mname1 = "JHI-Hv50k-2016-73692")
plotPXG(HvDRR48per, "JHI-Hv50k-2016-73692")

find.marker(HvDRR48per, 4, 41.7)
effectplot(HvDRR48per, mname1 = "JHI-Hv50k-2016-228890")
summary(HvDRR48per, mname1 = "JHI-Hv50k-2016-228890")
plotPXG(HvDRR48per, "JHI-Hv50k-2016-228890")

summary(fitqtl(HvDRR48per, qtl = HvDRR48rqtl, method = "hk", get.ests=TRUE, formula=Flowering~Q1 + Q2))

HvDRR48CI2 <- lodint(HvDRR48rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
HvDRR48CI4 <- lodint(HvDRR48rqtl, chr = 4, qtl.index = 2, expandtomarkers = TRUE)

pt(-4.990, 81, lower.tail = FALSE)
#p<0.01
pt(4.932, 81, lower.tail = FALSE)
#p<0.01
