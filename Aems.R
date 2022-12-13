dir()
source("main_function_aem_h2.R")
library(Matrix)
library(lme4)
library(emmeans)
setwd("/Aems/")

Phenotypes <- read.csv("/All_environments_FT.csv")

for(i in 1:nrow(Phenotypes)){
  for(j in 7:ncol(Phenotypes)){
    if(is.na(Phenotypes[i,j])){
      Phenotypes[i,j] <- 0.0
    }
  }
}



Phenotypes <- na.omit(Phenotypes)



Phenotypes$combine <- paste(Phenotypes$Genotype, Phenotypes$Environment, sep = "_")

uni.com <- unique(Phenotypes$combine)
length(uni.com)

table(Phenotypes$combine)[table(Phenotypes$combine) > 1]

select.aem <- NULL
for(i in 1:length(uni.com)){
  tmp <- Phenotypes[Phenotypes$combine == uni.com[i], ]
  if(length(unique(tmp$GenotypeAwnAverage)) == 1){
    select.aem <- rbind(select.aem, tmp[1,])
  }else{
    print(uni.com[i])
    stop("check dataset")
  }
}
rm(tmp)


traitnames <- colnames(Phenotypes)[5]

j <- 1

cat("trait:", traitnames[j], "\n")

formu.aem.origin <- as.formula(paste(traitnames[j], "~ (1|Environment) + Genotype"))
formu.h2.origin <- as.formula(paste(traitnames[j], " ~ (1|Environment) + (1|Genotype)"))

fit.aem.origin <- lmer(formu.aem.origin, data = select.awn)

fit.h2.origin <- lmer(formu.h2.origin, data = select.awn)
out.origin <- AEM.fun(fit = fit.aem.origin, select.name = "Genotype")
AEM.my.origin <- out.origin$aem
colnames(AEM.my.origin)[1] <- "Genotype"
write.csv(AEM.my.origin, file = 'AEMS_FT.csv')
errorvar.my.origin <- out.origin$contrast.var.mean/2

vcov.origin <- as.data.frame(VarCorr(fit.h2.origin))
vg.origin <- vcov.origin[vcov.origin$grp == "Genotype", "vcov"]

# heritability
h2 <- vg.origin/(vg.origin + errorvar.my.origin)
print(h2)
write.csv(h2, file = "h2_FT")

## fit model: G-fixed E-random G:E random

formu.aem.G_E <- as.formula(paste(traitnames[j], "~ (1|Environment) + Genotype + (1|Genotype:Environment)"))

fit.aem.G_E <- lmer(formu.aem.G_E, data = select.awn)
out.origin.G_E <- AEM.fun(fit = fit.aem.G_E, select.name = "Genotype")
AEM.my.origin_G_E <- out.origin.G_E$aem
colnames(AEM.my.origin_G_E)[1] <- "Genotype"
write.csv(AEM.my.origin_G_E, file = "AEM_G_E_FT")




dir()
source("main_function_aem_h2.R")
library(Matrix)
library(lme4)
library(emmeans)
setwd("/Aems/")

Phenotypes <- read.csv("/All_environments_PH.csv")

for(i in 1:nrow(Phenotypes)){
  for(j in 7:ncol(Phenotypes)){
    if(is.na(Phenotypes[i,j])){
      Phenotypes[i,j] <- 0.0
    }
  }
}



Phenotypes <- na.omit(Phenotypes)



Phenotypes$combine <- paste(Phenotypes$Genotype, Phenotypes$Environment, sep = "_")

uni.com <- unique(Phenotypes$combine)
length(uni.com)

table(Phenotypes$combine)[table(Phenotypes$combine) > 1]

select.aem <- NULL
for(i in 1:length(uni.com)){
  tmp <- Phenotypes[Phenotypes$combine == uni.com[i], ]
  if(length(unique(tmp$GenotypeAwnAverage)) == 1){
    select.aem <- rbind(select.aem, tmp[1,])
  }else{
    print(uni.com[i])
    stop("check dataset")
  }
}
rm(tmp)


traitnames <- colnames(Phenotypes)[5]

j <- 1

cat("trait:", traitnames[j], "\n")

formu.aem.origin <- as.formula(paste(traitnames[j], "~ (1|Environment) + Genotype"))
formu.h2.origin <- as.formula(paste(traitnames[j], " ~ (1|Environment) + (1|Genotype)"))

fit.aem.origin <- lmer(formu.aem.origin, data = select.awn)

fit.h2.origin <- lmer(formu.h2.origin, data = select.awn)
out.origin <- AEM.fun(fit = fit.aem.origin, select.name = "Genotype")
AEM.my.origin <- out.origin$aem
colnames(AEM.my.origin)[1] <- "Genotype"
write.csv(AEM.my.origin, file = 'AEMS_PH.csv')
errorvar.my.origin <- out.origin$contrast.var.mean/2

vcov.origin <- as.data.frame(VarCorr(fit.h2.origin))
vg.origin <- vcov.origin[vcov.origin$grp == "Genotype", "vcov"]

# heritability
h2 <- vg.origin/(vg.origin + errorvar.my.origin)
print(h2)
write.csv(h2, file = "h2_PH")

## fit model: G-fixed E-random G:E random

formu.aem.G_E <- as.formula(paste(traitnames[j], "~ (1|Environment) + Genotype + (1|Genotype:Environment)"))

fit.aem.G_E <- lmer(formu.aem.G_E, data = select.awn)
out.origin.G_E <- AEM.fun(fit = fit.aem.G_E, select.name = "Genotype")
AEM.my.origin_G_E <- out.origin.G_E$aem
colnames(AEM.my.origin_G_E)[1] <- "Genotype"
write.csv(AEM.my.origin_G_E, file = 'AEM_G_E_PH"
