#BMI

rm(list = ls())
library(readxl);
library(readr);
library(sva);
require('SmartSVA');
require('nlme');
library(parallel)
library(dplyr)

###CpG sites
load('E:/reviews/2024.1.11 HXM-BMI-heritability/data-mk/BMIz.Rdata')
bcpg <- bcpg[!duplicated(bcpg)]

###methylation data
load('E:/reviews/2024.1.11 HXM-BMI-heritability/data-mk/CorrectedBeta.Rdata')
###pheno data
PHENO <- read.csv("E:/reviews/2024.1.11 HXM-BMI-heritability/data-mk/pheno.csv",header=T,row.names = 1,as.is=T)

'%notin%' <- function(x,y)! (x %in% y)
which(bcpg%notin%rownames(beta))
sum(bcpg%in%rownames(beta))
bcpg <- bcpg[bcpg %in% rownames(beta)]
beta <- beta[bcpg,]

PHENO <- subset(PHENO,!is.na(PHENO$BMI))
PHENO <- subset(PHENO,!is.na(PHENO$sex))

PHENO_mz <- subset(PHENO,PHENO$zyg=='MZ')

myBeta.quantile_test <- beta
'%notin%' <- function(x,y)! (x %in% y)
PHENO <- PHENO[order(PHENO$Basename),]

sum(PHENO$Basename %notin% colnames(myBeta.quantile_test))
which(PHENO$Basename %notin% colnames(myBeta.quantile_test))
PHENO <- PHENO[-which(PHENO$Basename %notin% colnames(myBeta.quantile_test)),]
PHENOf <- data.frame(table(PHENO$TID))

table(PHENOf$Freq)
uniPHENO <- PHENO[PHENO$TID %in% PHENOf$Var[PHENOf$Freq != 2],]
PHENO <- PHENO[PHENO$TID %notin% uniPHENO$TID,]
myBeta.quantile_test <- myBeta.quantile_test[,which(colnames(myBeta.quantile_test)%in%PHENO$Basename)]
dim(PHENO)

na.myBeta.quantile_test <- na.omit(myBeta.quantile_test) 
dim(na.myBeta.quantile_test)
t_bn <- t(na.myBeta.quantile_test)
t_bn <- t_bn[order(rownames(t_bn)),] 
PHENO <- PHENO[which(PHENO$Basename %in% rownames(t_bn)),]
dim(PHENO)
rownames(PHENO) <- PHENO$Basename
PHENO <- PHENO[with(PHENO,order(rownames(t_bn))),]  
sum(rownames(PHENO)!=rownames(t_bn))   # 0

dim(t_bn) #1060， 2594
dim(PHENO) #1060，79

AGE <- PHENO[,"age"]
SEX <- PHENO[,"sex"]
smk <- PHENO[,"smk"]
drk <- PHENO[,"drk"]
BMI <- PHENO[,"BMI"]

fid <- as.factor(PHENO[,'TID'])
fid <- as.numeric(fid)

Y.r <- t(resid(lm(t_bn ~ BMI+AGE+SEX+smk+drk))) 
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1    
#（2）full model matrix
mod <- model.matrix(~ BMI+AGE+SEX+smk+drk) 
dim(mod)
svTable <- sva(t(t_bn), mod = mod, mod0 = NULL, n.sv = n.sv, method="irw")
sv <- as.matrix(svTable$sv)
dim(sv)  
rownames(sv) <- rownames(t_bn)
n = dim(sv)[2];n 
colnames(sv) <- paste("sv",1:n,sep='')

################## full population analysis #######################################
lmeFit <- function(Design, GeneSet, random = NULL, cl = detectCores()-1) {
  # Design = design; GeneSet = dasenB[1:10,];  random = c('fid');cl = detectCores(logical = F)
  if (!is.null(random)) {
    rdm <- Design[, random]
    Design <- Design[, colnames(Design) %notin% random]
    Zygn <- Design[,'ZygMZ']
    Design <- Design[,-which(colnames(Design)=='ZygMZ')]
  }
  
  Design <- Design[,-1]
  Xs <- paste(colnames(Design), collapse = '+')
  Formula <- as.formula(paste0('y ~', Xs))
  cat('The formula is:', sep = '\n')
  cat(as.character(Formula), sep = '')
  cat(' ', sep = '\n')
  
  Design <- as.data.frame(Design)
  paired <- rdm[duplicated(rdm)]
  Design <- cbind(Design,Zygn, rdm)
  Design <- Design[Design[,'rdm'] %in% paired, ]
  rdm <- rdm[rdm %in% paired]
  Zygn <- Zygn[rdm %in% paired]
  cat(paste('Test contans', nrow(Design), 'samples, selected samples are...'), seq = '\n')
  print(rownames(Design))
  GeneSet <- GeneSet[,rownames(Design)]
  
  u <- 1:nrow(GeneSet)
  cl <- makeCluster(cl)
  # parallel::clusterExport(cl,varlist = c("GeneSet", "Design", "Formula"))
  cat('Processing regression...', sep = '\n')
  
  res <- parSapply(
    cl,
    X = u,
    FUN = function(x) {
      if (!require(nlme))
        require(nlme, lib = "~/Rlib")
      i <- x
      temp <- Design
      temp$y <- as.numeric(GeneSet[i, ])
      temp <- temp[!is.na(temp$y),]
      index <- temp$rdm[duplicated(temp$rdm)]
      temp <- temp[temp$rdm %in% index,]
      tem <- rep(NA, 5)
      try({
        lm2 <- lme(Formula, random = ~1|rdm+1|Zygn, data = temp)
        tem <- summary(lm2)[['tTable']][2, ]
      }
      )
      return(tem)
    },
    simplify = T
  )
  stopCluster(cl)
  
  cat('The regression process is done', seq = '\n')
  
  res <- t(res)
  rownames(res) <- rownames(GeneSet)
  colnames(res) <- c("Value", "Std.Error", "DF", "t-value", "p-value")
  res <- res[order(res[, 'p-value']), ]
  res <- as.data.frame(res)
  res$p.adjust <- p.adjust(res$`p-value`, method = 'BH')
  
  result <- list()
  result$Design <- Design
  result$Formula <- Formula
  result$Coef <- res
  return(result)
}

#（6）Design Geneset
PHENO$ZygMZ <- ifelse(PHENO$Zyg=='MZ',0,1)

formula <- as.formula("~BMI+AGE+SEX+smk+drk+TID+ZygMZ")
design1 <- model.matrix(formula,PHENO) 
design1 <- cbind(design1,sv[rownames(design1),])

na.myBeta.quantile_test1 <- na.myBeta.quantile_test
design1 <- design1[rownames(design1) %in% colnames(na.myBeta.quantile_test1),]

#（7）mix model
result_BMI_full <- lmeFit(Design = design1, GeneSet = na.myBeta.quantile_test1, random = 'TID')
result_BMI_full_Coef <- result_BMI_full$Coef
sum(result_BMI_full_Coef$p.adjust<0.05,na.rm = T)


##################### MZ analysis #################################
lmeFit_cotwin <- function(Design, GeneSet, random = NULL, cl = detectCores()-1) {
  # Design = design; GeneSet = dasenB[1:10,];  random = c('fid');cl = detectCores(logical = F)
  if (!is.null(random)) {
    rdm <- Design[, random]
    Design <- Design[, colnames(Design) %notin% random]
  }
  
  Design <- Design[,-1]
  Xs <- paste(colnames(Design), collapse = '+')
  Formula <- as.formula(paste0('y ~', Xs))
  cat('The formula is:', sep = '\n')
  cat(as.character(Formula), sep = '')
  cat(' ', sep = '\n')
  
  Design <- as.data.frame(Design)
  paired <- rdm[duplicated(rdm)]
  Design <- cbind(Design, rdm)
  Design <- Design[Design[,'rdm'] %in% paired, ]
  rdm <- rdm[rdm %in% paired]
  cat(paste('Test contans', nrow(Design), 'samples, selected samples are...'), seq = '\n')
  print(rownames(Design))
  GeneSet <- GeneSet[,rownames(Design)]
  
  u <- 1:nrow(GeneSet)
  cl <- makeCluster(cl)
  # parallel::clusterExport(cl,varlist = c("GeneSet", "Design", "Formula"))
  cat('Processing regression...', sep = '\n')
  
  res <- parSapply(
    cl,
    X = u,
    FUN = function(x) {
      if (!require(nlme))
        require(nlme, lib = "~/Rlib")
      i <- x
      temp <- Design
      temp$y <- as.numeric(GeneSet[i, ])
      temp <- temp[!is.na(temp$y),]
      index <- temp$rdm[duplicated(temp$rdm)]
      temp <- temp[temp$rdm %in% index,]
      tem <- rep(NA, 5)
      try({
        lm2 <- lme(Formula, random = ~1|rdm, data = temp)
        tem <- summary(lm2)[['tTable']][2, ]
      }
      )
      return(tem)
    },
    simplify = T
  )
  stopCluster(cl)
  
  cat('The regression process is done', seq = '\n')
  
  res <- t(res)
  rownames(res) <- rownames(GeneSet)
  colnames(res) <- c("Value", "Std.Error", "DF", "t-value", "p-value")
  res <- res[order(res[, 'p-value']), ]
  res <- as.data.frame(res)
  res$p.adjust <- p.adjust(res$`p-value`, method = 'BH')
  
  result <- list()
  result$Design <- Design
  result$Formula <- Formula
  result$Coef <- res
  return(result)
}

###
PHENO_mz <- PHENO_mz[order(PHENO_mz$ID20),]
sum(PHENO_mz$Basename %notin% colnames(myBeta.quantile_test))
PHENO_mz <- PHENO_mz[-which(PHENO_mz$Basename %notin% colnames(myBeta.quantile_test)),]

PHENO_mzf <- data.frame(table(PHENO_mz$TID))
table(PHENO_mzf$Freq)   
PHENO_mzF <- PHENO_mzf[PHENO_mzf$Freq=="1",]
PHENO_mz <- PHENO_mz[PHENO_mz$TID %notin% PHENO_mzF$Var1,]  

#Beta
myBeta.quantile_test <- myBeta.quantile_test[,colnames(myBeta.quantile_test)%in%PHENO_mz$Basename]
na.myBeta.quantile_test <- na.omit(myBeta.quantile_test) ##delete 'NA' probes, save as sva analysis
dim(na.myBeta.quantile_test)  

t_beta <- t(na.myBeta.quantile_test)

PHENO_mz <- PHENO_mz[which(PHENO_mz$Basename %in% rownames(t_beta)),]
rownames(PHENO_mz) <- PHENO_mz$Basename
PHENO_mz <- PHENO_mz[with(PHENO_mz,order(ID20)),]   ##sort PHENO_mz DATA by Barcode

#make sure the same order of t_beta and na.myBeta.quantile_test with pheno_mz
t_beta <- t_beta[PHENO_mz$Basename,]
na.myBeta.quantile_test <- t(t_beta)

##sva adjust
platedata <- read.csv('E:/研究生阶段/复核/2024.1.11 HXM-BMI-heritability/data-mk/Plate.csv')
platedata <- subset(platedata,select=c('plate','Basename'))
rownames_1 <- rownames(t_beta)

#change the order of platedata 
rownames(platedata) <- platedata$Basename 
platedata <- platedata[rownames_1,]

#cbind t_beta and platedata
t_beta <- cbind(t_beta,rownames_1)
t_beta <- merge(t_beta,platedata,by.x='rownames_1',by.y='Basename', all.x = T)
batch_1 <- t_beta$plate
PHENO_mz <- PHENO_mz[which(PHENO_mz$Basename %in% t_beta$rownames_1),]

modcombat_1 <- model.matrix( ~ sex + age + smk + drk + BMI, data = PHENO_mz)
combat_data <- ComBat(dat = na.myBeta.quantile_test, batch = batch_1, mod = modcombat_1,par.prior = F,
                      prior.plots = FALSE,
                      mean.only = FALSE,
                      ref.batch = NULL,
                      BPPARAM = bpparam("SerialParam"))

beta_final <- t(combat_data)  
PHENO_mz <- PHENO_mz[order(PHENO_mz$ID20),]
beta_final <- beta_final[PHENO_mz$Basename,]

##twin difference
PHENO_mz$BMI <- as.numeric(PHENO_mz$BMI)
for(i in seq(1,length(PHENO_mz[,1]),2)){
  PHENO_mz[i,'BMI_mean'] <- mean(c(PHENO_mz[i,'BMI'],PHENO_mz[i+1,'BMI']),trim=0,na.rm=T)
  PHENO_mz[i+1,'BMI_mean'] <- mean(c(PHENO_mz[i,'BMI'],PHENO_mz[i+1,'BMI']),trim=0,na.rm=T)
  PHENO_mz[i,'BMI'] <- (PHENO_mz[i,'BMI']-mean(c(PHENO_mz[i,'BMI'],PHENO_mz[i+1,'BMI']),trim=0,na.rm=T))
  PHENO_mz[i+1,'BMI'] <- (-1)*PHENO_mz[i,'BMI']
}
mean_diff <- function(cpg,data){
  for (i in 1:length(cpg)) {
    cpg1 <- cpg[i]
    for(j in seq(1,length(data[,1]),2)){
      data[j,cpg1] <- (data[j,cpg1]-mean(c(data[j,cpg1],data[j+1,cpg1]),trim=0,na.rm=T))
      data[j+1,cpg1] <- (-1)*data[j,cpg1]
    }
  }
  return(data)
}


beta_final <- mean_diff(colnames(beta_final),beta_final)


fid <- as.factor(PHENO_mz[,'TID'])
fid <- as.numeric(fid)

smk <- PHENO_mz[,"smk"]
drk <- PHENO_mz[,"drk"]
BMI <- PHENO_mz[,"BMI"]

formula <- as.formula("~BMI+age+sex + smk+drk+fid")
design1 <- model.matrix(formula,PHENO_mz)
na.myBeta.quantile_test1 <- t(beta_final)
design1 <- design1[rownames(design1) %in% colnames(na.myBeta.quantile_test1),]

#get the result
#MZ
result_BMI_mz <- lmeFit_cotwin(Design = design1, GeneSet = na.myBeta.quantile_test1, random = 'fid')
result_BMI_mz_Coef <- result_BMI_mz$Coef
result_BMI_mz <- subset(result_BMI_mz_Coef,result_BMI_mz_Coef$p.adjust<0.05)
dim(result_BMI_mz)
cpg_BMI_mz <- rownames(result_BMI_mz)

#full population
result_BMI_full_Coef <- subset(result_BMI_full_Coef,result_BMI_full_Coef$p.adjust<0.05)
dim(result_BMI_full_Coef);
dim(result_BMI_mz)
cpg_BMI_f <- rownames(result_BMI_full_Coef)

#CpG
cpg_BMI <- c(cpg_BMI_f,cpg_BMI_mz)
cpg_BMI <- cpg_BMI[!duplicated(cpg_BMI)]


####
load('E:/revies/2024.1.11 HXM-BMI-heritability/data-mk/Annotation.RData')
Annotation <- subset(Annotation,select=c(Name,CHR,MAPINFO,UCSC_RefGene_Name,UCSC_RefGene_Group,Enhancer))

result_BMI_full_Coef$CpG <- rownames(result_BMI_full_Coef)
BMI_res <- merge(result_BMI_full_Coef,Annotation,by.x='CpG',by.y='Name',all.x=T)

#change the meaning of variance
BMI_res$UCSC_RefGene_Name <- ifelse(BMI_res$UCSC_RefGene_Name=="",'unannotated',BMI_res$UCSC_RefGene_Name)
BMI_res$UCSC_RefGene_Group <- ifelse(BMI_res$UCSC_RefGene_Group=="",'unannotated',BMI_res$UCSC_RefGene_Group)
BMI_res$dir <- ifelse(BMI_res$Value <0 ,'-','+')
BMI_res <- subset(BMI_res,select=c(CpG,CHR,MAPINFO,UCSC_RefGene_Name,UCSC_RefGene_Group,Enhancer,dir,Value,p.adjust))
rownames(BMI_res) <- BMI_res$CpG
BMI_res <- BMI_res[rownames(result_BMI_full_Coef),]
cpg_BMI <- rownames(result_BMI_full_Coef)



load('E:/reviews/2024.1.11 HXM-BMI-heritability/data-mk/BMIz.Rdata')


bcpgz_summary <- bcpgz %>%
  group_by(cpg) %>%
  summarise(
    non_na_count = sum(!is.na(beta))
  )

bcpgz <- bcpgz %>%
  left_join(bcpgz_summary, by = "cpg") %>%
  mutate(
    beta = as.character(case_when(
      non_na_count == 0 ~ "NA",
      non_na_count == 1 ~ ifelse(beta > 0, "+", "-"),
      any(beta > 0) && any(beta < 0) ~ "+-",
      all(beta > 0) ~ "+",
      all(beta < 0) ~ "-",
      TRUE ~ "NA"
    ))
  ) %>%
  select(-non_na_count)

bcpgz <- bcpgz %>% distinct(cpg, .keep_all = TRUE)

BMI_res <- merge(BMI_res,bcpgz,by.x = 'CpG',by.y = 'cpg',all.x = T)

number <- which((BMI_res$dir=='-'&!is.na(BMI_res$beta)&BMI_res$beta=='+')|(BMI_res$dir=='+'&!is.na(BMI_res$beta)&BMI_res$beta=='-'))
table(number)



#delete discordancy
cpg_BMI <- cpg_BMI[-number] 
BMI_final <- subset(BMI_res, BMI_res$CpG %in% cpg_BMI)
      
      
      ######## annotation
      load('E:/reviews/2024.1.11 HXM-BMI-heritability/data-mk/Annotation.RData')
      Annotation <- subset(Annotation,select=c(Name,CHR,MAPINFO,UCSC_RefGene_Name,UCSC_RefGene_Group,Enhancer))
      result_BMI_mz$CpG <- rownames(result_BMI_mz)
      BMI_mz_res <- merge(result_BMI_mz,Annotation,by.x='CpG',by.y='Name',all.x=T)
      
      BMI_mz_res$UCSC_RefGene_Name <- ifelse(BMI_mz_res$UCSC_RefGene_Name=="",'unannotated',BMI_mz_res$UCSC_RefGene_Name)
      BMI_mz_res$UCSC_RefGene_Group <- ifelse(BMI_mz_res$UCSC_RefGene_Group=="",'unannotated',BMI_mz_res$UCSC_RefGene_Group)
      BMI_mz_res$dir <- ifelse(BMI_mz_res$Value<0,'-','+')
      
      BMI_mz_res <- subset(BMI_mz_res,select=c(CpG,CHR,MAPINFO,UCSC_RefGene_Name,UCSC_RefGene_Group,Enhancer,dir,Value,p.adjust))
      rownames(BMI_mz_res) <- BMI_mz_res$CpG
      BMI_mz_res <- BMI_mz_res[rownames(result_BMI_mz),]
      cpg_BMI_mz <- rownames(result_BMI_mz)
      
      
      library(dplyr)
      bcpgz_summary <- bcpgz %>%
        group_by(cpg) %>%
        summarise(
          non_na_count = sum(!is.na(beta))
        )
      
      bcpgz <- bcpgz %>%
        left_join(bcpgz_summary, by = "cpg") %>%
        mutate(
          beta = as.character(case_when(
            non_na_count == 0 ~ "NA",
            non_na_count == 1 ~ ifelse(beta > 0, "+", "-"),
            any(beta > 0) && any(beta < 0) ~ "+-",
            all(beta > 0) ~ "+",
            all(beta < 0) ~ "-",
            TRUE ~ "NA"
          ))
        ) %>%
        select(-non_na_count)
      
      bcpgz <- bcpgz %>% distinct(cpg, .keep_all = TRUE)
      BMI_mz_res <- merge(BMI_mz_res,bcpgz,by.x = 'CpG',by.y = 'cpg',all.x = T)
   
number_mz <- which((BMI_mz_res$dir=='-'&!is.na(BMI_mz_res$beta)&BMI_mz_res$beta=='+')|(BMI_mz_res$dir=='+'&!is.na(BMI_mz_res$beta)&BMI_mz_res$beta=='-'))

number_mz

      
      
      
      
      
      
      