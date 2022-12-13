#FLOWERING TIME####

setwd("/QTLanalysis/MPP_analysis/")


#GENOTYPIC DATA####
SNP50K <- read.csv("MPP_analysis/A.3_clean_data_across_populations.csv")
SNP50K <- data.frame(lapply(SNP50K, as.character), stringsAsFactors=FALSE)

colnames(SNP50K) <- gsub("X", "", colnames(SNP50K), fixed = TRUE)
colnames(SNP50K) <- paste(substr(colnames(SNP50K), 1,5), "-", substr(colnames(SNP50K), 6,7), "-", 
                       substr(colnames(SNP50K), 8,12), "-", substr(colnames(SNP50K), 13,16), sep= "")
rownames(SNP50K) <- SNP50K[,1]
SNP50K <- SNP50K[,-1]
SNP50K <- as.matrix(SNP50K)

source("libraryfieldexperimentation.R")
Genotypes <- read.csv("Genotypes.csv", sep = ";")
Genotypes <- data.frame(lapply(Genotypes, as.character), stringsAsFactors=FALSE)

for (i in 1:ncol(SNP50K)){
  colnames(SNP50K)[i] <- Genotypes$Genotypecode[which(Genotypes$Lotcode == colnames(SNP50K)[i])]
}

SNP50K[SNP50K!= "GG" & SNP50K != "AA" & SNP50K != "CC" & SNP50K != "TT"] <- NA

SNP50K <- t(SNP50K)
SNP50K <- SNP50K[!duplicated(rownames(SNP50K)), ]

write.csv(SNP50K, file = "MPP_Geno.csv")
MPP_Geno <- read.csv("MPP_Geno.csv")

MPP_Geno <- data.frame(lapply(MPP_Geno, as.character), stringsAsFactors=FALSE)

colnames(MPP_Geno) <- MPP_Geno[1,]
MPP_Geno <- MPP_Geno[-1,]

rownames(MPP_Geno) <- MPP_Geno[,1]
MPP_Geno <- MPP_Geno[,-1]

MPP_Geno <- SNP50K
MPP_Geno <- data.frame(lapply(MPP_Geno, as.character), stringsAsFactors=FALSE)
MPP_Geno <- MPP_Geno[!duplicated(MPP_Geno$X), ]
rownames(MPP_Geno) <- MPP_Geno[,1]
MPP_Geno <- MPP_Geno[,-1]

Parents <- read.csv("B.1_parent_color_table.csv")
Parents <- data.frame(lapply(Parents, as.character), stringsAsFactors=FALSE)
Parents <- as.data.frame(Parents$parent)
Parents <- data.frame(lapply(Parents, as.character), stringsAsFactors=FALSE)

geno.off <- SNP50K[-which(row.names(SNP50K) %in% unique(Parents$Parents.parent)),]
geno.par <- SNP50K[which(row.names(SNP50K) %in% unique(Parents$Parents.parent)),]

write.csv(SNP50K, file = "MPP_Geno.csv")
write.csv(geno.off, file = "Geno_off.csv")
write.csv(geno.par, file = "Geno_par.csv")

Geno_off <- read.csv("Geno_off.csv", row.names = 1)
Geno_par <- read.csv("Geno_par.csv", row.names = 1)

#####

#MAP DATA####
Chr1cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr1.txt", sep = " ", stringsAsFactors = FALSE)
Chr2cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr2.txt", sep = " ", stringsAsFactors = FALSE)
Chr3cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr3.txt", sep = " ", stringsAsFactors = FALSE)
Chr4cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr4.txt", sep = " ", stringsAsFactors = FALSE)
Chr5cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr5.txt", sep = " ", stringsAsFactors = FALSE)
Chr6cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr6.txt", sep = " ", stringsAsFactors = FALSE)
Chr7cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr7.txt", sep = " ", stringsAsFactors = FALSE)

Consensus_map <- rbind(Chr1cons,Chr2cons,Chr3cons,Chr4cons,Chr5cons,Chr6cons,Chr7cons)
Consensus_map <- Consensus_map[,-3]
Consensus_map <- na.omit(Consensus_map)
rownames(Consensus_map) <- 1:nrow(Consensus_map)
Consensus_map <- data.frame(lapply(Consensus_map, as.character), stringsAsFactors=FALSE)
Consensus_map$chr <- as.integer(as.numeric(substr(Consensus_map$chr,4,4)))
Consensus_map$cM <- as.numeric(Consensus_map$cM)
colnames(Geno_off) <- gsub(".", "-", colnames(Geno_off), fixed = TRUE)
colnames(Geno_par) <- gsub(".", "-", colnames(Geno_par), fixed = TRUE)

Consensus_map <- Consensus_map[which(Consensus_map[,1] %in% colnames(Geno_off)),]
Consensus_map <- Consensus_map[which(Consensus_map[,1] %in% colnames(Geno_par)),]

rownames(Consensus_map) <- 1:nrow(Consensus_map)
colnames(Consensus_map) <- c("mk.names", "chr", "pos.cM")

par_clu <- read.csv("Genetic_haplotypes.csv", row.names = 1)
Consensus_map <- Consensus_map[which(Consensus_map$mk.names %in% rownames(par_clu)),]

Geno_par <- Geno_par[,which(colnames(Geno_par) %in% Consensus_map[,1])]
Geno_off <- Geno_off[,which(colnames(Geno_off) %in% Consensus_map[,1])]

set_rownames_par <- rownames(Geno_par)
set_rownames_off <- rownames(Geno_off)

Geno_par <- data.frame(lapply(Geno_par, as.character), stringsAsFactors=FALSE)
Geno_par <- as.matrix(Geno_par)
rownames(Geno_par) <- set_rownames_par


Geno_off <- data.frame(lapply(Geno_off, as.character), stringsAsFactors=FALSE)
Geno_off <- as.matrix(Geno_off)
rownames(Geno_off) <- set_rownames_off

#####

#PHENOTYPIC DATA####

aems <- read.csv("AEM_origin_alltraits_allpop_17to19_FT_PH_LA.csv")
aems <- data.frame(lapply(aems, as.character), stringsAsFactors=FALSE)

source("libraryfieldexperimentation.R")
Genotypes <- read.csv("Genotypes.csv", sep = ";")
Genotypes <- data.frame(lapply(Genotypes, as.character), stringsAsFactors=FALSE)

Genotype_Barley_AM_Bkilian <- Genotypes[which(Genotypes$Population == "Barley_AM_Bkilian"),]
Geno_off <- Geno_off[-which(rownames(Geno_off) %in% Genotype_Barley_AM_Bkilian$Genotypecode),]
set_rownames_off <- rownames(Geno_off)


for (i in 1:nrow(Geno_off)){
  Lotcode <- as.data.frame(Genotypes$Lotcode[which(Genotypes$Genotypecode == rownames(Geno_off)[i])])
  Lotcode <- as.character(Lotcode[1,])
  rownames(Geno_off)[i]<- gobacktoF3(Lotcode, Genotypes)
}

Pheno <- aems[which(aems$ind %in% rownames(Geno_off)),]

matches_for_aems_ind <- as.data.frame(cbind(rownames(Geno_off), set_rownames_off))
matches_for_aems_ind <- data.frame(lapply(matches_for_aems_ind, as.character), stringsAsFactors=FALSE)

for (i in 1:nrow(Pheno)){
  Pheno$ind[i] <- matches_for_aems_ind$set_rownames_off[which(matches_for_aems_ind$V1 == Pheno$ind[i])]
}

row.names(Geno_off) <- set_rownames_off


rownames(Pheno) <- Pheno$ind
colnames(Pheno)[4] <- "Flowering"
Pheno <- Pheno[,4, drop=FALSE]

Pheno$Flowering <- as.numeric(Pheno$Flowering)
Pheno <- as.matrix(Pheno)

#####

#CROSS PARENT INFORMATION####

Crosses <- read.csv("C.2_parents_A_B_info.csv", row.names = 1)
Crosses <- data.frame(lapply(Crosses, as.character), stringsAsFactors=FALSE)
Crosses <- as.matrix(Crosses)
attr(Crosses, "dimnames") <- NULL

#DATA PROCESSING####

colnames(Geno_off) <- gsub(".", "-", colnames(Geno_off), fixed = TRUE)
colnames(Geno_par) <- gsub(".", "-", colnames(Geno_par), fixed = TRUE)

markers_vector <- Consensus_map[,1]
Geno_par <- Geno_par[,markers_vector]
Geno_off <- Geno_off[,markers_vector]

#to order the columns of Geno_par and Geno_off with the same order of the 1st column of Consensus map (that contains the markers)

pheno_vector <- row.names(Geno_off)[which(rownames(Geno_off) %in% rownames(Pheno))]

Pheno <- Pheno[order(match(rownames(Pheno), pheno_vector)), , drop = FALSE]
Geno_off <- Geno_off[which(rownames(Geno_off) %in% rownames(Pheno)),]

#to order the rows of Pheno with the same order of the rows of Geno_off and then select in Geno_off only the genotypes for
#which I have phenotypic information

#####

#CROSS INDIVIDUALS####
Cross.ind <- rownames(Geno_off)

Genotypes <- read.csv("Genotypes.csv", sep = ";")
Genotypes$Population <- as.character(Genotypes$Population)

for (j in 1:length(Cross.ind)){
  info_pop <- as.data.frame(Genotypes$Population[which(Genotypes$Genotypecode == Cross.ind[j])])
  Cross.ind[j] <- as.character(info_pop[1,])
}

print(unique(Cross.ind))
length(Cross.ind)

#####

#DATA PROCESSING####
library(mppR)

print(unique(Cross.ind))
print(Crosses)

nrow(Geno_off)
length(Cross.ind)

mppData <- create.mppData(geno.off = Geno_off, geno.par = Geno_par, map = Consensus_map, pheno = Pheno, 
                          cross.ind = Cross.ind, par.per.cross = Crosses)

saveRDS(mppData, file = "mppData_ancestral_with_GHap_FT_genetic_map.RDS")



library(mppR)
library(ggplot2)
library(farver)
my.loc <- ("/QTLanalysis/MPP_analysis/")
mppData <- readRDS("mppData_ancestral_with_GHap_FT_genetic_map.RDS")
summary(mppData)
attr(mppData$geno.off, "dimnames")[[2]]

par_clu <- read.csv("Genetic_haplotypes.csv", row.names = 1)
colnames(par_clu)[8] <- "W23829/803911"
colnames(par_clu)[19] <- "Unumli-Arpa"
par_clu <- as.matrix(par_clu)


mppData <- QC.mppData(mppData = mppData, n.lim = 15, MAF.pop.lim = 0.05, MAF.cr.miss = TRUE, mk.miss = 0.1,
                      gen.miss = 0.25, verbose = TRUE)

mppData <- IBS.mppData(mppData = mppData)

mppData <- IBD.mppData(mppData = mppData, het.miss.par = TRUE, type = "RIL", type.mating = "selfing")

#par_clu <- read.csv("/gpfs/project/projects/qggp/Francesco_barley/1stpaper/Rprojects/statistics/QTLanalysis/MPP_analysis/QTL_effect_ancestral/par_clu.csv", row.names = 1)
#colnames(par_clu)[8] <- "W23829/803911"
#colnames(par_clu)[19] <- "Unumli-Arpa"
#par_clu <- as.matrix(par_clu)

test <- attr(mppData$geno.off, "dimnames"[1])[[2]]
par_clu <- par_clu[which(rownames(par_clu) %in% test),]
#It is necessary to execute this two lines because after the quality control (QC.mppData) some markers from the mppData are deleted
#and it happens that par.clu and mppData don't match anymore, for this reason it is first necessary to delete those markers also
#from par.clu

mppdata_vector <- as.data.frame(attr(mppData$geno.off, "dimnames")[[2]], stringsAsFactors = FALSE)
class(attr(mppData$geno.off, "dimnames")[[2]])
par_clu_vector <- as.data.frame(row.names(par_clu), stringsAsFactors = FALSE)
class(rownames(par_clu))

check <- cbind(mppdata_vector, par_clu_vector)
check[3] <- ""
colnames(check)[3] <- "control"

for (i in 1:nrow(check)){
  if(as.character(check[i,1]) == check[i,2]){
    check[i,3] <- "YES"
  } else {
    check[i,3] <- "NO"
  }
}

print(check)




mppData <- parent_cluster.mppData(mppData = mppData, par.clu = par_clu)

summary(mppData)

SIM <- mpp_SIM(mppData = mppData, Q.eff = "anc")

SIM

cofactors <- QTL_select(Qprof = SIM)

cofactors

CIM <- mpp_CIM(mppData = mppData, Q.eff = "anc", cofactors = cofactors, plot.gen.eff = TRUE)

CIM

QTL <- QTL_select(Qprof = CIM)

QTL

gen.eff <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "anc")

gen.eff

summary(gen.eff, QTL = 1)

plot(x = CIM, QTL = QTL, type = "l")

Plot_result <- plot(x = CIM, QTL = QTL, type = "l")

save(Plot_result, file = "Plot_result_ancestral_GHap_FT_genetic_map.RData")

plot(x = CIM, gen.eff= TRUE, mppData = mppData, QTL = QTL, Q.eff = "anc", main = "QTL genetic effect plot")

Plot_parents_effect <- plot(x = CIM, gen.eff= TRUE, mppData = mppData, QTL = QTL, Q.eff = "anc", main = "QTL genetic effect plot")

save(Plot_parents_effect, file = "Plot_parents_effect_ancestral_GHap_FT_genetic_map.RData")


CV <- mpp_CV(pop.name = "HvDRR", trait.name = "Flowering", mppData = mppData, Q.eff = "anc", her = 0.4, Rep = 1,
             k = 3, verbose = FALSE, output.loc = my.loc)

save(CV, file = "CV_ancestral_GHap_FT_genetic_map.RData")



#PLANT HEIGHT####

setwd("/QTLanalysis/MPP_analysis/")

#GENOTYPIC DATA####
SNP50K <- read.csv("A.3_clean_data_across_populations.csv")
SNP50K <- data.frame(lapply(SNP50K, as.character), stringsAsFactors=FALSE)

colnames(SNP50K) <- gsub("X", "", colnames(SNP50K), fixed = TRUE)
colnames(SNP50K) <- paste(substr(colnames(SNP50K), 1,5), "-", substr(colnames(SNP50K), 6,7), "-", 
                       substr(colnames(SNP50K), 8,12), "-", substr(colnames(SNP50K), 13,16), sep= "")
rownames(SNP50K) <- SNP50K[,1]
SNP50K <- SNP50K[,-1]
NP50K <- as.matrix(SNP50K)

source("libraryfieldexperimentation.R")
Genotypes <- read.csv("Genotypes.csv", sep = ";")
Genotypes <- data.frame(lapply(Genotypes, as.character), stringsAsFactors=FALSE)

for (i in 1:ncol(SNP50K)){
  colnames(SNP50K)[i] <- Genotypes$Genotypecode[which(Genotypes$Lotcode == colnames(SNP50K)[i])]
}

SNP50K[SNP50K!= "GG" & SNP50K != "AA" & SNP50K != "CC" & SNP50K != "TT"] <- NA

SNP50K <- t(SNP50K)
SNP50K <- SNP50K[!duplicated(rownames(SNP50K)), ]

write.csv(SNP50K, file = "MPP_Geno.csv")
MPP_Geno <- read.csv("MPP_Geno.csv")

MPP_Geno <- data.frame(lapply(MPP_Geno, as.character), stringsAsFactors=FALSE)

colnames(MPP_Geno) <- MPP_Geno[1,]
MPP_Geno <- MPP_Geno[-1,]

rownames(MPP_Geno) <- MPP_Geno[,1]
MPP_Geno <- MPP_Geno[,-1]

MPP_Geno <- SNP50K
MPP_Geno <- data.frame(lapply(MPP_Geno, as.character), stringsAsFactors=FALSE)
MPP_Geno <- MPP_Geno[!duplicated(MPP_Geno$X), ]
rownames(MPP_Geno) <- MPP_Geno[,1]
MPP_Geno <- MPP_Geno[,-1]

Parents <- read.csv("B.1_parent_color_table.csv")
Parents <- data.frame(lapply(Parents, as.character), stringsAsFactors=FALSE)
Parents <- as.data.frame(Parents$parent)
Parents <- data.frame(lapply(Parents, as.character), stringsAsFactors=FALSE)

geno.off <- SNP50K[-which(row.names(SNP50K) %in% unique(Parents$Parents.parent)),]
geno.par <- SNP50K[which(row.names(SNP50K) %in% unique(Parents$Parents.parent)),]

write.csv(SNP50K, file = "MPP_Geno.csv")
write.csv(geno.off, file = "MPP_analysis/Geno_off.csv")
write.csv(geno.par, file = "MPP_analysis/Geno_par.csv")

Geno_off <- read.csv("MPP_analysis/Geno_off.csv", row.names = 1)
Geno_par <- read.csv("MPP_analysis/Geno_par.csv", row.names = 1)

#####

#MAP DATA####
Chr1cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr1.txt", sep = " ", stringsAsFactors = FALSE)
Chr2cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr2.txt", sep = " ", stringsAsFactors = FALSE)
Chr3cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr3.txt", sep = " ", stringsAsFactors = FALSE)
Chr4cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr4.txt", sep = " ", stringsAsFactors = FALSE)
Chr5cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr5.txt", sep = " ", stringsAsFactors = FALSE)
Chr6cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr6.txt", sep = " ", stringsAsFactors = FALSE)
Chr7cons <- read.csv("Integreted_Map_based_on_DRR_RecRate_DVI/Integrated_map_chr7.txt", sep = " ", stringsAsFactors = FALSE)

Consensus_map <- rbind(Chr1cons,Chr2cons,Chr3cons,Chr4cons,Chr5cons,Chr6cons,Chr7cons)
Consensus_map <- Consensus_map[,-3]
Consensus_map <- na.omit(Consensus_map)
rownames(Consensus_map) <- 1:nrow(Consensus_map)
Consensus_map <- data.frame(lapply(Consensus_map, as.character), stringsAsFactors=FALSE)
Consensus_map$chr <- as.integer(as.numeric(substr(Consensus_map$chr,4,4)))
Consensus_map$cM <- as.numeric(Consensus_map$cM)
colnames(Geno_off) <- gsub(".", "-", colnames(Geno_off), fixed = TRUE)
colnames(Geno_par) <- gsub(".", "-", colnames(Geno_par), fixed = TRUE)

Consensus_map <- Consensus_map[which(Consensus_map[,1] %in% colnames(Geno_off)),]
Consensus_map <- Consensus_map[which(Consensus_map[,1] %in% colnames(Geno_par)),]

rownames(Consensus_map) <- 1:nrow(Consensus_map)
colnames(Consensus_map) <- c("mk.names", "chr", "pos.cM")
#to set up the map with the data requested by the mppR package

par_clu <- read.csv("Genetic_haplotypes.csv", row.names = 1)
Consensus_map <- Consensus_map[which(Consensus_map$mk.names %in% rownames(par_clu)),]
#to exclude the markers that are not in par_clu

Geno_par <- Geno_par[,which(colnames(Geno_par) %in% Consensus_map[,1])]
Geno_off <- Geno_off[,which(colnames(Geno_off) %in% Consensus_map[,1])]

#Here I initially selected only the markers of Consensus_map that are in included in the Geno files and then the other way around
#in order to have same markers indentifiers for the different datasets

set_rownames_par <- rownames(Geno_par)
set_rownames_off <- rownames(Geno_off)

Geno_par <- data.frame(lapply(Geno_par, as.character), stringsAsFactors=FALSE)
Geno_par <- as.matrix(Geno_par)
rownames(Geno_par) <- set_rownames_par


Geno_off <- data.frame(lapply(Geno_off, as.character), stringsAsFactors=FALSE)
Geno_off <- as.matrix(Geno_off)
rownames(Geno_off) <- set_rownames_off

#####

#PHENOTYPIC DATA####

aems <- read.csv("AEM_origin_alltraits_allpop_17to19_FT_PH_LA.csv")
aems <- data.frame(lapply(aems, as.character), stringsAsFactors=FALSE)

source("libraryfieldexperimentation.R")
Genotypes <- read.csv("Genotypes.csv", sep = ";")
Genotypes <- data.frame(lapply(Genotypes, as.character), stringsAsFactors=FALSE)

Genotype_Barley_AM_Bkilian <- Genotypes[which(Genotypes$Population == "Barley_AM_Bkilian"),]
Geno_off <- Geno_off[-which(rownames(Geno_off) %in% Genotype_Barley_AM_Bkilian$Genotypecode),]
set_rownames_off <- rownames(Geno_off)
#to delete the genotypes that don't belong to any of the 45 populations

for (i in 1:nrow(Geno_off)){
  Lotcode <- as.data.frame(Genotypes$Lotcode[which(Genotypes$Genotypecode == rownames(Geno_off)[i])])
  Lotcode <- as.character(Lotcode[1,])
  rownames(Geno_off)[i]<- gobacktoF3(Lotcode, Genotypes)
}

Pheno <- aems[which(aems$ind %in% rownames(Geno_off)),]

matches_for_aems_ind <- as.data.frame(cbind(rownames(Geno_off), set_rownames_off))
matches_for_aems_ind <- data.frame(lapply(matches_for_aems_ind, as.character), stringsAsFactors=FALSE)

for (i in 1:nrow(Pheno)){
  Pheno$ind[i] <- matches_for_aems_ind$set_rownames_off[which(matches_for_aems_ind$V1 == Pheno$ind[i])]
}

row.names(Geno_off) <- set_rownames_off


rownames(Pheno) <- Pheno$ind
colnames(Pheno)[3] <- "Plant_Height"
Pheno <- Pheno[,3, drop=FALSE]

Pheno$Plant_Height <- as.numeric(Pheno$Plant_Height)
Pheno <- as.matrix(Pheno)

#####

#CROSS PARENT INFORMATION####

Crosses <- read.csv("C.2_parents_A_B_info.csv", row.names = 1)
Crosses <- data.frame(lapply(Crosses, as.character), stringsAsFactors=FALSE)
Crosses <- as.matrix(Crosses)
attr(Crosses, "dimnames") <- NULL

#DATA PROCESSING####

colnames(Geno_off) <- gsub(".", "-", colnames(Geno_off), fixed = TRUE)
colnames(Geno_par) <- gsub(".", "-", colnames(Geno_par), fixed = TRUE)

markers_vector <- Consensus_map[,1]
Geno_par <- Geno_par[,markers_vector]
Geno_off <- Geno_off[,markers_vector]

#to order the columns of Geno_par and Geno_off with the same order of the 1st column of Consensus map (that contains the markers)

pheno_vector <- row.names(Geno_off)[which(rownames(Geno_off) %in% rownames(Pheno))]

Pheno <- Pheno[order(match(rownames(Pheno), pheno_vector)), , drop = FALSE]
Geno_off <- Geno_off[which(rownames(Geno_off) %in% rownames(Pheno)),]

#to order the rows of Pheno with the same order of the rows of Geno_off and then select in Geno_off only the genotypes for
#which I have phenotypic information

#####

#CROSS INDIVIDUALS####
Cross.ind <- rownames(Geno_off)

Genotypes <- read.csv("Genotypes.csv", sep = ";")
Genotypes$Population <- as.character(Genotypes$Population)

for (j in 1:length(Cross.ind)){
  info_pop <- as.data.frame(Genotypes$Population[which(Genotypes$Genotypecode == Cross.ind[j])])
  Cross.ind[j] <- as.character(info_pop[1,])
}

print(unique(Cross.ind))
length(Cross.ind)

#####

#DATA PROCESSING####
library(mppR)

print(unique(Cross.ind))
print(Crosses)

nrow(Geno_off)
length(Cross.ind)

mppData <- create.mppData(geno.off = Geno_off, geno.par = Geno_par, map = Consensus_map, pheno = Pheno, 
                          cross.ind = Cross.ind, par.per.cross = Crosses)

saveRDS(mppData, file = "mppData_ancestral_with_GHap_PH_genetic_map.RDS")


library(mppR)
library(ggplot2)
library(farver)
my.loc <- ("QTLanalysis/MPP_analysis/")
mppData <- readRDS("mppData_ancestral_with_GHap_PH_genetic_map.RDS")
summary(mppData)
attr(mppData$geno.off, "dimnames")[[2]]

par_clu <- read.csv("Genetic_haplotypes.csv", row.names = 1)
colnames(par_clu)[8] <- "W23829/803911"
colnames(par_clu)[19] <- "Unumli-Arpa"
par_clu <- as.matrix(par_clu)


mppData <- QC.mppData(mppData = mppData, n.lim = 15, MAF.pop.lim = 0.05, MAF.cr.miss = TRUE, mk.miss = 0.1,
                      gen.miss = 0.25, verbose = TRUE)

mppData <- IBS.mppData(mppData = mppData)

mppData <- IBD.mppData(mppData = mppData, het.miss.par = TRUE, type = "RIL", type.mating = "selfing")

par_clu <- read.csv("QTL_effect_ancestral/par_clu.csv", row.names = 1)
colnames(par_clu)[8] <- "W23829/803911"
colnames(par_clu)[19] <- "Unumli-Arpa"
par_clu <- as.matrix(par_clu)

test <- attr(mppData$geno.off, "dimnames"[1])[[2]]
par_clu <- par_clu[which(rownames(par_clu) %in% test),]

mppdata_vector <- as.data.frame(attr(mppData$geno.off, "dimnames")[[2]], stringsAsFactors = FALSE)
class(attr(mppData$geno.off, "dimnames")[[2]])
par_clu_vector <- as.data.frame(row.names(par_clu), stringsAsFactors = FALSE)
class(rownames(par_clu))

check <- cbind(mppdata_vector, par_clu_vector)
check[3] <- ""
colnames(check)[3] <- "control"

for (i in 1:nrow(check)){
  if(as.character(check[i,1]) == check[i,2]){
    check[i,3] <- "YES"
  } else {
    check[i,3] <- "NO"
  }
}

print(check)




mppData <- parent_cluster.mppData(mppData = mppData, par.clu = par_clu)

summary(mppData)

SIM <- mpp_SIM(mppData = mppData, Q.eff = "anc")

SIM

cofactors <- QTL_select(Qprof = SIM)

cofactors

CIM <- mpp_CIM(mppData = mppData, Q.eff = "anc", cofactors = cofactors, plot.gen.eff = TRUE)

CIM

QTL <- QTL_select(Qprof = CIM)

QTL

gen.eff <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "anc")

gen.eff

summary(gen.eff, QTL = 1)

plot(x = CIM, QTL = QTL, type = "l")

Plot_result <- plot(x = CIM, QTL = QTL, type = "l")

save(Plot_result, file = "Plot_result_ancestral_GHap_PH_genetic_map.RData")

plot(x = CIM, gen.eff= TRUE, mppData = mppData, QTL = QTL, Q.eff = "anc", main = "QTL genetic effect plot")

Plot_parents_effect <- plot(x = CIM, gen.eff= TRUE, mppData = mppData, QTL = QTL, Q.eff = "anc", main = "QTL genetic effect plot")

save(Plot_parents_effect, file = "Plot_parents_effect_ancestral_GHap_PH_genetic_map.RData")


CV <- mpp_CV(pop.name = "HvDRR", trait.name = "Flowering", mppData = mppData, Q.eff = "anc", her = 0.4, Rep = 1,
             k = 3, verbose = FALSE, output.loc = my.loc)

save(CV, file = "CV_ancestral_GHap_PH_genetic_map.RData")



