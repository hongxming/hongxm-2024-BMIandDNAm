rm(list = ls())
library(OpenMx,lib='~/Rlib')
mxOption( NULL, 'Default optimizer' ,'SLSQP' )
mxOption( NULL, "mvnRelEps", value= mxOption(NULL, key="mvnRelEps")/5)
getOption('mxOptions')$mvnRelEps
#load('/public/home/hongxming/epigenoheri/deta.Rdata')
load('/public/home/hongxming/data_longi/deta_longi2.Rdata')
PHENO_2 <- subset(PHENO_2,select=c(ID20,BMI_1))
deta <- merge(deta,PHENO_2,by.x='ID20_t1',by.y='ID20',all.x=T)
deta <- merge(deta,PHENO_2,by.x='ID20_t2',by.y='ID20',all.x=T)
deta$BMI_t1 <- deta$BMI_1.x
deta$BMI_t2 <- deta$BMI_1.y
deta_re <- deta
load('/public/home/hongxming/bmi_ace/result_frame.Rdata')
ACESummer_re <- ACESummer
AESummer_re <- AESummer
dropRcSummer_re <- dropRcSummer
result_bi$A11lo <- NA
result_bi$A11 <- NA
result_bi$A11up <- NA
result_bi$A12lo <- NA
result_bi$A12 <- NA
result_bi$A12up <- NA
result_bi$A22lo <- NA
result_bi$A22 <- NA
result_bi$A22up <- NA
result_bi$C11lo <- NA
result_bi$C11 <- NA
result_bi$C11up <- NA
result_bi$C12lo <- NA
result_bi$C12 <- NA
result_bi$C12up <- NA
result_bi$C22lo <- NA
result_bi$C22 <- NA
result_bi$C22up <- NA
result_bi$E11lo <- NA
result_bi$E11 <- NA
result_bi$E11up <- NA
result_bi$E12lo <- NA
result_bi$E12 <- NA
result_bi$E12up <- NA
result_bi$E22lo <- NA
result_bi$E22 <- NA
result_bi$E22up <- NA
result_bi$ACEAIC <- NA
result <- result_bi

# PHENO_bmi <- subset(PHENO,select = c(ID20,BMI))
# deta <- merge(deta,PHENO_bmi,by.x='ID20_t1',by.y='ID20',all.x = T)
# deta <- merge(deta,PHENO_bmi,by.x='ID20_t2',by.y='ID20',all.x = T)

load('/public/home/hongxming/bmi_ace/cpg_bmi_all.Rdata')
cpg_bmi_all <- cpg_bmi_all[which(cpg_bmi_all%in%CPG_all)]
CPG_all <- cpg_bmi_all

# deta <- subset(deta,!is.na(deta$bw))
###########################################################
###########################################################
#set variables' parameters

for (i in 1:length(CPG_all)) {
    tryCatch({
    deta <- deta_re
    vart <- CPG_all[i]
    
    ###########################################################
    ###########################################################
    #set variables' parameters
    vars      <- c('BMI', vart)         # list of variables names
    nv        <- 2                               # number of variables
    ntv       <- nv*2                            # number of total variables
    covars    <- c('age_2', 'sex')
    ncv       <- length(covars)
    
    selVars   <- paste(vars,   c(rep(1,nv), rep(2,nv)), sep = "_t")
    covVars   <- paste(covars, c(rep(1,nv), rep(2,nv)), sep = "_t")
    ageVars <- c("age_2_t1","age_2_t2")
    
    deta[,ageVars] = scale(deta[,ageVars])
    mzData  <- subset(deta, Zyg_t1=='MZ', c(selVars,covVars))
    dzData  <- subset(deta, Zyg_t1=='DZ', c(selVars,covVars))
    sum(is.na(mzData))
    sum(is.na(dzData))
    mzData <- na.omit(mzData)
    dzData <- na.omit(dzData)
    dim(mzData); dim(dzData)
    mzData    <- mzData[complete.cases(mzData), c(selVars,covVars)]
    dzData    <- dzData[complete.cases(dzData), c(selVars,covVars)]
    
    
    #------------ace/threshold crude estimate---------------------------
    r_MZ <- cor(mzData[,1:2],use='complete.obs')
    r_MZ
    r_DZ <- cor(dzData[,1:2],use='complete.obs')
    r_DZ
    
    ### Mean
    meanGMZ <- c(mean(unlist(mzData[,c(1,3)])), mean(unlist(mzData[,c(2,4)])))
    meanGDZ <- c(mean(unlist(dzData[,c(1,3)])), mean(unlist(dzData[,c(2,4)])))
    
    ### Covariance Matrices
    CovMaMZ <- matrix(NA, ntv, ntv)
    CovMZ <- c()
    for(i in 1:ntv){
      for(j in i:ntv){
        CovMaMZ[i,j] <- cov(mzData[,i], mzData[j])
        CovMaMZ[j,i] <- CovMaMZ[i,j]
        CovMZ <- c(CovMZ, CovMaMZ[i,j])
      }
    }
    CovMaDZ <- matrix(NA, ntv, ntv)
    CovDZ <- c()
    for(i in 1:ntv){
      for(j in i:ntv){
        CovMaDZ[i,j] <- cov(dzData[,i], dzData[j])
        CovMaDZ[j,i] <- CovMaDZ[i,j]
        CovDZ <- c(CovDZ, CovMaDZ[i,j])
      }
    }
    lbvas <- matrix(NA, ntv, ntv)
    diag(lbvas) <- min(min(CovMaMZ), min(CovMaDZ))
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Create Labels
    labCovar <- covVars[c(1,3,2,4)]
    labBeta  <- cbind(paste("B", covars, "_", vars[1], sep = ""), paste("B", covars, "_", vars[2], sep = ""))
    
    laMeMZT <- c("BMI1MZT", "HBA1C1MZT", "BMI2MZT", "HBA1C2MZT")                   
    laMeDZT <- c("BMI1DZT", "HBA1C1DZT", "BMI2DZT", "HBA1C2DZT")
    
    laVaMZ    <- c("MZVP1T1", "MZCVT1", "MZWP1", "MZCTT1T2","MZVP2T1", "MZCTT2T1", "MZWP2", "MZVP1T2", "MZCVT2", "MZVP2T2")  # labels for (co)variances for MZ twins
    laVaDZ    <- c("DZVP1T1", "DZCVT1", "DZWP1", "DZCTT1T2","DZVP2T1", "DZCTT2T1", "DZWP2", "DZVP1T2", "DZCVT2", "DZVP2T2")  # labels for (co)variances for DZ twins
    
    # ------------------------------------------------------------------------------
    #                            Twin 1                           Twin 2
    #                  Phenotype 1  |  Phenotype 2  ||  Phenotype 1  |  Phenotype 2
    #        
    # Twin 1 Pheno 1       MZVP1T1
    #        Pheno 2        MZCVT1         MZVP2T1
    #        
    # Twin 2 Pheno 1         MZWP1        MZCTT2T1          MZVP1T2
    #        Pheno 2      MZCTT1T2           MZWP2           MZCVT2         MZVP2T2 
    # ------------------------------------------------------------------------------
    
    # ----------------------------------------------------------------------------------------------------------------------
    # PREPARE MODEL
    # Create Algebra for expected Mean Matrices
    covarsMZ   <- mxMatrix(type="Full", nrow=2, ncol=ncv, free=F, label=labCovar, name="MZDefVars")
    covarsDZ   <- mxMatrix(type="Full", nrow=2, ncol=ncv, free=F, label=labCovar, name="DZDefVars")
    covbeta    <- mxMatrix(type="Full", nrow=ncv, ncol=nv, free=T, values=.1, labels=labBeta, name="Beta")
    
    meanMZ   <- mxMatrix(type="Full", nrow=1, ncol=ntv, free=TRUE, values=meanGMZ, labels=laMeMZT, name="meanMZ" )
    meanDZ   <- mxMatrix(type="Full", nrow=1, ncol=ntv, free=TRUE, values=meanGDZ, labels=laMeDZT, name="meanDZ" )
    expMeanMZ <- mxAlgebra(expression = meanMZ+ cbind(MZDefVars[1,]%*%Beta, MZDefVars[2,]%*%Beta), name="expMeanMZ")
    expMeanDZ <- mxAlgebra(expression = meanDZ+ cbind(DZDefVars[1,]%*%Beta, DZDefVars[2,]%*%Beta), name="expMeanDZ")
    
    # Create Algebra for expected Variance/Covariance Matrices
    covMZ <- mxMatrix(type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=CovMZ, lbound=lbvas, labels=laVaMZ, name="covMZ" )
    covDZ <- mxMatrix(type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=CovDZ, lbound=lbvas, labels=laVaDZ, name="covDZ" )
    
    ### 为啥要有个相关阵？
    corMZ <- mxAlgebra(expression = cov2cor(CovMaMZ), name = "corMZ")
    corDZ <- mxAlgebra(expression = cov2cor(CovMaDZ), name = "corDZ")
    
    # Create Data Objects for Multiple Groups
    dataMZ <- mxData(observed = mzData, type = "raw")
    dataDZ <- mxData(observed = dzData, type = "raw")
    
    # Create Expectation Objects for Multiple Groups
    expMZ <- mxExpectationNormal(covariance="covMZ", means="meanMZ", dimnames=selVars)
    expDZ <- mxExpectationNormal(covariance="covDZ", means="meanDZ", dimnames=selVars)
    funML <- mxFitFunctionML()
    
    # Create Model Objects for Multiple Groups
    modelMZ <- mxModel(covarsMZ, covbeta, meanMZ, expMeanMZ, covMZ, corMZ, dataMZ, expMZ, funML, name="MZ")
    modelDZ <- mxModel(covarsDZ, covbeta, meanDZ, expMeanDZ, covDZ, corDZ, dataDZ, expDZ, funML, name="DZ")
    multi <- mxFitFunctionMultigroup(c("MZ","DZ"))
    
    # Create Confidence Interval Objects
    ciCov <- mxCI(c('MZ.covMZ','DZ.covDZ'))
    ciMean <- mxCI(c('MZ.meanMZ','DZ.meanDZ'))
    ciCor <- mxCI(c("MZ.corMZ", "DZ.corDZ"))
    
    # Build Saturated Model with Confidence Intervals
    modelSAT <- mxModel( "twoSATc", modelMZ, modelDZ, multi, ciCov, ciMean, ciCor)
    
    # ----------------------------------------------------------------------------------------------------------------------
    # RUN MODEL
    # Run Saturated Model
    fitSAT <- mxTryHard(modelSAT, intervals = FALSE, extraTries = 15)
    sumSAT <- summary(fitSAT)
    

    # if(r_MZ[1,2]<2*r_DZ[1,2]){
    svPa <- .2 # start value for path coefficient
    svPe <- .5 # start value for path coefficient for e
    svRa <- .2 # start value for genetic correlations
    svRc <- .2 # start value for shared env correlations
    svRe <- .1 # start value for unique env correlations
    
    # ----------------------------------------------------------------------------------------------------------------------
    # PREPARE MODEL
    # Create Algebra for expected Mean Matrices
    labLower <- function(lab,nv)    {paste(lab,rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="")}
    labVars  <- function(lab,vars)  {paste(lab,vars,sep="")}
    labSdiag <- function(lab,nv)    {paste(lab,rev(nv+1-sequence(1:(nv-1))),rep(1:(nv-1),(nv-1):1),sep="") }
    labOdiag <- function(lab,nv)    {paste(lab,c(rev(nv+1-sequence(1:(nv-1))),rep(1:(nv-1),(nv-1):1)),c(rep(1:(nv-1),(nv-1):1),rev(nv+1-sequence(1:(nv-1)))),sep="")}
    labDiag  <- function(lab,nv)    {paste(lab,1:nv,1:nv,sep="")} 
    labFull  <- function(lab,nr,nc) {paste(lab,1:nr,rep(1:nc,each=nr),sep="")}
    
    covars   <- mxMatrix(type="Full", nrow=2, ncol=ncv, free=F, label=labCovar, name="DefVars")
    meanZ    <- mxMatrix(type="Full", nrow=1, ncol=ntv, free=TRUE, values=meanGMZ, labels=labVars("mean",vars), name="meanZ")
    expMeanZ <- mxAlgebra(expression= meanZ + cbind(DefVars[1,]%*%Beta, DefVars[1,]%*%Beta), name="expMeanZ") #+ cbind(bm%*%Age,bm%*%Age)
    
    # Create Matrices for Path Coefficients
    pathA <- mxMatrix(type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svPa, label=labDiag("a",nv), lbound=.0001, name="a")
    pathC <- mxMatrix(type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svPa, label=labDiag("c",nv), lbound=.0001, name="c")
    pathE <- mxMatrix(type="Diag", nrow=nv, ncol=nv, free=TRUE, values=svPe, label=labDiag("e",nv), lbound=.0001, name="e")
    
    # Create Matrices for Correlation Coefficients within/across Individuals
    vart_t1 <- paste0(vart,'_t1')
    vart_t2 <- paste0(vart,'_t2')
    if((cor(deta[,vart_t1],deta[,'BMI_t1'])+cor(deta[,vart_t2],deta[,'BMI_t2']))>=0){
    pathRa <- mxMatrix(type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRa, label=labSdiag("ra",nv), lbound=0, ubound=1, name="Ra")
    pathRc <- mxMatrix(type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRc, label=labSdiag("rc",nv), lbound=0, ubound=1, name="Rc")
    pathRe <- mxMatrix(type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRe, label=labSdiag("re",nv), lbound=0, ubound=1, name="Re")
    }else{
    pathRa <- mxMatrix(type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRa, label=labSdiag("ra",nv), lbound=-1, ubound=0, name="Ra")
    pathRc <- mxMatrix(type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRc, label=labSdiag("rc",nv), lbound=-1, ubound=0, name="Rc")
    pathRe <- mxMatrix(type="Stand", nrow=nv, ncol=nv, free=TRUE, values=svRe, label=labSdiag("re",nv), lbound=-1, ubound=0, name="Re")
    }
    # Create Matrices for Causal Path Coefficients - Causal Direction: down column, out row
    ### 这是个啥呀...
    pathB <- mxMatrix(type="Full", nrow=nv, ncol=nv, free=c(F,F,F,F), labels=labFull("b",nv,nv), lbound=-.99, ubound=.99, name="b")
    
    # # Algebra to Constrain Correlation Matrices to be Positive Definite
    # cnstrPos <- mxMatrix(type="Full", nrow=1, ncol=6, free=FALSE, values=.0001, name="cnstrPos")
    # corMin <- mxAlgebra(expression= cbind(min(eigenval(Ra)), min(eigenval(Rc)), min(eigenval(Re))), name="corMin")
    # corPos <- mxConstraint(expression= corMin > cnstrPos, name="corPos")
    
    # Create Algebra for Variance Components
    covA <- mxAlgebra(expression= a %*% (Ra) %*% t(a), name="A")
    covC <- mxAlgebra(expression= c %*% (Rc) %*% t(c), name="C")
    covE <- mxAlgebra(expression= e %*% (Re) %*% t(e), name="E")
    
    # Create Algebra for Causal Variance Components
    matI <- mxMatrix(type="Iden", nrow=nv, ncol=nv, name="I")
    invIB <- mxAlgebra(expression= solve(I-b), name="IB")
    
    # Create Algebra for Total Variances and Standard Deviations (diagonal only)
    covV  <- mxAlgebra(expression= A+C+E, name="V")
    invSD <- mxAlgebra(expression= solve(sqrt(I*V)), name="iSD")
    covIV <- mxAlgebra(expression= IB %&% V, name="IV")
    
    # Create Algebra for Expected Variance/Covariance Matrices in MZ & DZ twins
    covMZ <- mxAlgebra(expression= A + C, name="cMZ")
    covDZ <- mxAlgebra(expression= 0.5%x%(A)+ C, name="cDZ")
    matI2 <- mxMatrix( type="Iden", nrow=2, ncol=2, name="I2")
    expCovMZ <- mxAlgebra( expression= I2 %x% IB %&% rbind(cbind(V, cMZ), cbind(cMZ, V)), name="expCovMZ")
    expCovDZ <- mxAlgebra( expression= I2 %x% IB %&% rbind(cbind(V, cDZ), cbind(cDZ, V)), name="expCovDZ")
    matZ2 <- mxMatrix( type="Zero", nrow=2, ncol=2, name="Z2")
    
    # Create Data Objects for Multiple Groups
    dataMZ <- mxData(observed=mzData, type="raw")
    dataDZ <- mxData(observed=dzData, type="raw")
    
    # Create Expectation Objects for Multiple Groups
    expMZ <- mxExpectationNormal(covariance="expCovMZ", means="expMeanZ", dimnames=selVars)
    expDZ <- mxExpectationNormal(covariance="expCovDZ", means="expMeanZ", dimnames=selVars)
    funML <- mxFitFunctionML()
    
    Rph <- mxAlgebra(expression = cov2cor(V), name = "Rph")
    ren <- mxAlgebra(expression = cov2cor(C+E), name = "ren")
    a2 <- mxAlgebra(expression = A/V, name = "a2")
    c2 <- mxAlgebra(expression = C/V, name = "c2")
    e2 <- mxAlgebra(expression = E/V, name = "e2")
    pa <- mxAlgebra(expression=sqrt(a2[1,1])*sqrt(a2[2,2])*Ra[1,2]/Rph[1,2],name = "pa")
    pc <- mxAlgebra(expression=sqrt(c2[1,1])*sqrt(c2[2,2])*Rc[1,2]/Rph[1,2],name = "pc")
    pe <- mxAlgebra(expression=sqrt(e2[1,1])*sqrt(e2[2,2])*Re[1,2]/Rph[1,2],name = "pe")
    
    ciACE <- mxCI(c("a[1,1]","a[2,2]","c[1,1]","c[2,2]","e[1,1]","e[2,2]","A[1,1]","A[2,1]","A[2,2]","C[1,1]","C[2,1]","C[2,2]",
                    "E[1,1]","E[2,1]","E[2,2]","a2[1,1]","a2[1,2]","a2[2,2]","c2[1,1]","c2[1,2]","c2[2,2]",
                    "e2[1,1]","e2[1,2]","e2[2,2]","Rph[1,2]","Ra[2,1]","Rc[2,1]","Re[2,1]","pa","pc","pe","ren[1,2]"))
    calc <- list(Rph, a2, c2, e2, pa, pc, pe, ren, ciACE)
    pars <- list(matI, matI2, matZ2)
    parsT <- list(meanZ)
    parsZ <- list(pathA, pathC, pathE, pathRa, pathRc, pathRe, covA, covC, covE, covV, invSD, pathB, invIB, covIV)
    modelMZ <- mxModel(pars, parsZ, covars, covbeta, meanZ, expMeanZ, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ")
    modelDZ <- mxModel(pars, parsZ, covars, covbeta, meanZ, expMeanZ, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ")
    modelT <- list(modelMZ, modelDZ)
    multi <- mxFitFunctionMultigroup(c("MZ","DZ"))
    modelACE <- mxModel("ACE", pars, parsT, parsZ, modelT, multi, calc)
    fitACE <- mxTryHard(modelACE, extraTries=15, greenOK=TRUE, silent=TRUE, intervals = TRUE)
    mxCompare(fitSAT, c(fitACE))  ### 0.3961476
    summary(fitACE)   ### AIC:991.951   BIC:1037.005

    
    ACESumm <- summary(fitACE,verbose = T)
    ACESumm$CI[16:32, ]; ACESumm$CIdetail[31:64,1:3]

    
    # ----------------------------------------------------------------------------------------------------------------------
    # RUN SUBMODELS
    
    # Run AE model
    dropC <- omxSetParameters(fitACE,labels=c("c11","c22"),free=F,values=0,name="dropC")
    # dropCFit  <- mxTryHard(dropC, intervals=T, extraTries = 15)
    # mxCompare(fitACE, c(dropCFit))  # yes #
    # summary(dropCFit,verbose = T)
    # dropCSumm <- summary(dropCFit,verbose = T)
    # dropCSumm$CI
    # dropCSumm$CIdetail[,1:3]
    
    AEfit <- omxSetParameters(dropC,labels=c("rc21"), free=F, values=0, name="AEfit")
    AEfit  <- mxTryHard(AEfit, intervals=T, extraTries = 15)
    # mxCompare(fitACE, c(dropCFit, AEfit))

    dropC1 <- omxSetParameters(fitACE,labels=c("c11"),free=F,values=0,name="dropC1")
    dropC2 <- omxSetParameters(fitACE,labels=c("c22"),free=F,values=0,name="dropC2")
    
    dropRcfit <- omxSetParameters(dropC1,labels=c("rc21"),free=F,values=0,name="dropRc")
    dropRcfit  <- mxTryHard(dropRcfit, intervals=T, extraTries = 15)
    dropRcSumm <- summary(dropRcfit,verbose = T)
    
    fitcompare <- mxCompare(fitACE, c(AEfit,dropRcfit))  # yes #
    mxCompare(fitSAT, c(fitACE))
    
    AESumm <- summary(AEfit,verbose = T)
    AESumm$CI[25:32, ]
    AESumm$CIdetail[49:64, 1:3]
    dropRcSumm <- summary(dropRcfit,verbose = T)
    dropRcSumm$CI[25:32, ]
    dropRcSumm$CIdetail[31:64, 1:3]
    
    ACESummer <- ACESummer_re
    AESummer <- AESummer_re
    dropRcSummer <- dropRcSummer_re
    
    if(is.null(ACESumm$CI[25,2])){
    fitACE <- mxTryHard(modelACE, extraTries=15, greenOK=TRUE, silent=TRUE, intervals = TRUE)
    ACESumm <- summary(fitACE,verbose = T)
    }
    if(is.null(AESumm$CI[25,2])){
    AEfit  <- mxTryHard(AEfit, intervals=T, extraTries = 15)
    AESumm <- summary(AEfit,verbose = T)
    }
    if(is.null(dropRcSumm$CI[25,2])){
    dropRcfit  <- mxTryHard(dropRcfit, intervals=T, extraTries = 15)
    dropRcSumm <- summary(dropRcfit,verbose = T)
    }

    if(is.null(ACESumm$CI[25,2])){
      fitACE <- mxTryHard(modelACE, extraTries=15, greenOK=TRUE, silent=TRUE, intervals = TRUE)
      ACESumm <- summary(fitACE,verbose = T)
    }
    if(is.null(AESumm$CI[25,2])){
      AEfit  <- mxTryHard(AEfit, intervals=T, extraTries = 15)
      AESumm <- summary(AEfit,verbose = T)
    }
    if(is.null(dropRcSumm$CI[25,2])){
      dropRcfit  <- mxTryHard(dropRcfit, intervals=T, extraTries = 15)
      dropRcSumm <- summary(dropRcfit,verbose = T)
    }
    
    if(!is.null(ACESumm$CI[25,2])){
      ACESumm$CI[7,1] <- ACESumm$CIdetail[13,3]
      ACESumm$CI[7,3] <- ACESumm$CIdetail[14,3]
      ACESumm$CI[8,1] <- ACESumm$CIdetail[15,3]
      ACESumm$CI[8,3] <- ACESumm$CIdetail[16,3]
      ACESumm$CI[9,1] <- ACESumm$CIdetail[17,3]
      ACESumm$CI[9,3] <- ACESumm$CIdetail[18,3]
      ACESumm$CI[10,1] <- ACESumm$CIdetail[19,3]
      ACESumm$CI[10,3] <- ACESumm$CIdetail[20,3]
      ACESumm$CI[11,1] <- ACESumm$CIdetail[21,3]
      ACESumm$CI[11,3] <- ACESumm$CIdetail[22,3]
      ACESumm$CI[12,1] <- ACESumm$CIdetail[23,3]
      ACESumm$CI[12,3] <- ACESumm$CIdetail[24,3]
      ACESumm$CI[13,1] <- ACESumm$CIdetail[25,3]
      ACESumm$CI[13,3] <- ACESumm$CIdetail[26,3]
      ACESumm$CI[14,1] <- ACESumm$CIdetail[27,3]
      ACESumm$CI[14,3] <- ACESumm$CIdetail[28,3]
      ACESumm$CI[15,1] <- ACESumm$CIdetail[29,3]
      ACESumm$CI[15,3] <- ACESumm$CIdetail[30,3]
      ACESumm$CI[25,1] <- ACESumm$CIdetail[49,3]
      ACESumm$CI[25,3] <- ACESumm$CIdetail[50,3]
      ACESumm$CI[26,1] <- ACESumm$CIdetail[51,3]
      ACESumm$CI[26,3] <- ACESumm$CIdetail[52,3]
      ACESumm$CI[27,1] <- ACESumm$CIdetail[53,3]
      ACESumm$CI[27,3] <- ACESumm$CIdetail[54,3]
      ACESumm$CI[28,1] <- ACESumm$CIdetail[55,3]
      ACESumm$CI[28,3] <- ACESumm$CIdetail[56,3]
      ACESumm$CI[29,1] <- ACESumm$CIdetail[57,3]
      ACESumm$CI[29,3] <- ACESumm$CIdetail[58,3]
      ACESumm$CI[30,1] <- ACESumm$CIdetail[59,3]
      ACESumm$CI[30,3] <- ACESumm$CIdetail[60,3]
      ACESumm$CI[31,1] <- ACESumm$CIdetail[61,3]
      ACESumm$CI[31,3] <- ACESumm$CIdetail[62,3]
      ACESumm$CI[32,1] <- ACESumm$CIdetail[63,3]
      ACESumm$CI[32,3] <- ACESumm$CIdetail[64,3]
      ACESummer[c(25:32,7:15),] <- ACESumm$CI[c(25:32,7:15),]
    }else{
      ACESummer <- ACESummer_re
    }
    if(!is.null(AESumm$CI[25,2])){
      AESumm$CI[7,1] <- AESumm$CIdetail[13,3]
      AESumm$CI[7,3] <- AESumm$CIdetail[14,3]
      AESumm$CI[8,1] <- AESumm$CIdetail[15,3]
      AESumm$CI[8,3] <- AESumm$CIdetail[16,3]
      AESumm$CI[9,1] <- AESumm$CIdetail[17,3]
      AESumm$CI[9,3] <- AESumm$CIdetail[18,3]
      AESumm$CI[10,1] <- AESumm$CIdetail[19,3]
      AESumm$CI[10,3] <- AESumm$CIdetail[20,3]
      AESumm$CI[11,1] <- AESumm$CIdetail[21,3]
      AESumm$CI[11,3] <- AESumm$CIdetail[22,3]
      AESumm$CI[12,1] <- AESumm$CIdetail[23,3]
      AESumm$CI[12,3] <- AESumm$CIdetail[24,3]
      AESumm$CI[13,1] <- AESumm$CIdetail[25,3]
      AESumm$CI[13,3] <- AESumm$CIdetail[26,3]
      AESumm$CI[14,1] <- AESumm$CIdetail[27,3]
      AESumm$CI[14,3] <- AESumm$CIdetail[28,3]
      AESumm$CI[15,1] <- AESumm$CIdetail[29,3]
      AESumm$CI[15,3] <- AESumm$CIdetail[30,3]
      AESumm$CI[25,1] <- AESumm$CIdetail[49,3]
      AESumm$CI[25,3] <- AESumm$CIdetail[50,3]
      AESumm$CI[26,1] <- AESumm$CIdetail[51,3]
      AESumm$CI[26,3] <- AESumm$CIdetail[52,3]
      AESumm$CI[27,1] <- AESumm$CIdetail[53,3]
      AESumm$CI[27,3] <- AESumm$CIdetail[54,3]
      AESumm$CI[28,1] <- AESumm$CIdetail[55,3]
      AESumm$CI[28,3] <- AESumm$CIdetail[56,3]
      AESumm$CI[29,1] <- AESumm$CIdetail[57,3]
      AESumm$CI[29,3] <- AESumm$CIdetail[58,3]
      AESumm$CI[30,1] <- AESumm$CIdetail[59,3]
      AESumm$CI[30,3] <- AESumm$CIdetail[60,3]
      AESumm$CI[31,1] <- AESumm$CIdetail[61,3]
      AESumm$CI[31,3] <- AESumm$CIdetail[62,3]
      AESumm$CI[32,1] <- AESumm$CIdetail[63,3]
      AESumm$CI[32,3] <- AESumm$CIdetail[64,3]
      AESummer[c(25:32,7:15),] <- AESumm$CI[c(25:32,7:15),]
    }else{
      AESummer <- AESummer_re
    }
    if(!is.null(dropRcSumm$CI[25,2])){
      dropRcSumm$CI[7,1] <- dropRcSumm$CIdetail[13,3]
      dropRcSumm$CI[7,3] <- dropRcSumm$CIdetail[14,3]
      dropRcSumm$CI[8,1] <- dropRcSumm$CIdetail[15,3]
      dropRcSumm$CI[8,3] <- dropRcSumm$CIdetail[16,3]
      dropRcSumm$CI[9,1] <- dropRcSumm$CIdetail[17,3]
      dropRcSumm$CI[9,3] <- dropRcSumm$CIdetail[18,3]
      dropRcSumm$CI[10,1] <- dropRcSumm$CIdetail[19,3]
      dropRcSumm$CI[10,3] <- dropRcSumm$CIdetail[20,3]
      dropRcSumm$CI[11,1] <- dropRcSumm$CIdetail[21,3]
      dropRcSumm$CI[11,3] <- dropRcSumm$CIdetail[22,3]
      dropRcSumm$CI[12,1] <- dropRcSumm$CIdetail[23,3]
      dropRcSumm$CI[12,3] <- dropRcSumm$CIdetail[24,3]
      dropRcSumm$CI[13,1] <- dropRcSumm$CIdetail[25,3]
      dropRcSumm$CI[13,3] <- dropRcSumm$CIdetail[26,3]
      dropRcSumm$CI[14,1] <- dropRcSumm$CIdetail[27,3]
      dropRcSumm$CI[14,3] <- dropRcSumm$CIdetail[28,3]
      dropRcSumm$CI[15,1] <- dropRcSumm$CIdetail[29,3]
      dropRcSumm$CI[15,3] <- dropRcSumm$CIdetail[30,3]
      dropRcSumm$CI[25,1] <- dropRcSumm$CIdetail[49,3]
      dropRcSumm$CI[25,3] <- dropRcSumm$CIdetail[50,3]
      dropRcSumm$CI[26,1] <- dropRcSumm$CIdetail[51,3]
      dropRcSumm$CI[26,3] <- dropRcSumm$CIdetail[52,3]
      dropRcSumm$CI[27,1] <- dropRcSumm$CIdetail[53,3]
      dropRcSumm$CI[27,3] <- dropRcSumm$CIdetail[54,3]
      dropRcSumm$CI[28,1] <- dropRcSumm$CIdetail[55,3]
      dropRcSumm$CI[28,3] <- dropRcSumm$CIdetail[56,3]
      dropRcSumm$CI[29,1] <- dropRcSumm$CIdetail[57,3]
      dropRcSumm$CI[29,3] <- dropRcSumm$CIdetail[58,3]
      dropRcSumm$CI[30,1] <- dropRcSumm$CIdetail[59,3]
      dropRcSumm$CI[30,3] <- dropRcSumm$CIdetail[60,3]
      dropRcSumm$CI[31,1] <- dropRcSumm$CIdetail[61,3]
      dropRcSumm$CI[31,3] <- dropRcSumm$CIdetail[62,3]
      dropRcSumm$CI[32,1] <- dropRcSumm$CIdetail[63,3]
      dropRcSumm$CI[32,3] <- dropRcSumm$CIdetail[64,3]
      dropRcSummer[c(25:32,7:15),] <- dropRcSumm$CI[c(25:32,7:15),]
    }else{
      dropRcSummer <- dropRcSummer_re
    }
      result[vart,'model1'] <- 'ACE'
      result[vart,'model2'] <- fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]
      result[vart,'ACEAIC']<- fitcompare[1,'AIC']
      ACESummer <- ACESummer[,-4]
      AESummer <- AESummer[,-4]
      dropRcSummer <- dropRcSummer[,-4]
#      ACESummer <- ACESummer[-4,]
#      AESummer <- AESummer[-4,]
#      dropRcSummer <- dropRcSummer[-4,]
      result[vart,4] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[25,1],
                               ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[25,1],round(dropRcSummer, 2 )[25,1]))
      result[vart,5] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[25,2],
                               ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[25,2],round(dropRcSummer, 2 )[25,2]))
      result[vart,6] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[25,3],
                               ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[25,3],round(dropRcSummer, 2 )[25,3]))
      result[vart,7] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[26,1],
                               ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[26,1],round(dropRcSummer, 2 )[26,1]))
      result[vart,8] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[26,2],
                               ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[26,2],round(dropRcSummer, 2 )[26,2]))
      result[vart,9] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[26,3],
                               ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[26,3],round(dropRcSummer, 2 )[26,3]))
      result[vart,10] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[27,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[27,1],round(dropRcSummer, 2 )[27,1]))
      result[vart,11] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[27,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[27,2],round(dropRcSummer, 2 )[27,2]))
      result[vart,12] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[27,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[27,3],round(dropRcSummer, 2 )[27,3]))
      result[vart,13] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[32,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[32,1],round(dropRcSummer, 2 )[32,1]))
      result[vart,14] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[32,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[32,2],round(dropRcSummer, 2 )[32,2]))
      result[vart,15] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[32,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[32,3],round(dropRcSummer, 2 )[32,3]))
      result[vart,16] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[29,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[29,1],round(dropRcSummer, 2 )[29,1]))
      result[vart,17] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[29,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[29,2],round(dropRcSummer, 2 )[29,2]))
      result[vart,18] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[29,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[29,3],round(dropRcSummer, 2 )[29,3]))
      result[vart,19] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[30,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[30,1],round(dropRcSummer, 2 )[30,1]))
      result[vart,20] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[30,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[30,2],round(dropRcSummer, 2 )[30,2]))
      result[vart,21] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[30,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[30,3],round(dropRcSummer, 2 )[30,3]))
      result[vart,22] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[31,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[31,1],round(dropRcSummer, 2 )[31,1]))
      result[vart,23] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[31,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[31,2],round(dropRcSummer, 2 )[31,2]))
      result[vart,24] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[31,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[31,3],round(dropRcSummer, 2 )[31,3]))
      result[vart,25] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[7,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[7,1],round(dropRcSummer, 2 )[7,1]))
      result[vart,26] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[7,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[7,2],round(dropRcSummer, 2 )[7,2]))
      result[vart,27] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[7,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[7,3],round(dropRcSummer, 2 )[7,3]))
      result[vart,28] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[8,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[8,1],round(dropRcSummer, 2 )[8,1]))
      result[vart,29] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[8,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[8,2],round(dropRcSummer, 2 )[8,2]))
      result[vart,30] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[8,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[8,3],round(dropRcSummer, 2 )[8,3]))
      result[vart,31] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[9,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[9,1],round(dropRcSummer, 2 )[9,1]))
      result[vart,32] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[9,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[9,2],round(dropRcSummer, 2 )[9,2]))
      result[vart,33] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[9,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[9,3],round(dropRcSummer, 2 )[9,3]))
      result[vart,34] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[10,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[10,1],round(dropRcSummer, 2 )[10,1]))
      result[vart,35] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[10,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[10,2],round(dropRcSummer, 2 )[10,2]))
      result[vart,36] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[10,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[10,3],round(dropRcSummer, 2 )[10,3]))
      result[vart,37] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[11,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[11,1],round(dropRcSummer, 2 )[11,1]))
      result[vart,38] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[11,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[11,2],round(dropRcSummer, 2 )[11,2]))
      result[vart,39] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[11,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[11,3],round(dropRcSummer, 2 )[11,3]))
      result[vart,40] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[12,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[12,1],round(dropRcSummer, 2 )[12,1]))
      result[vart,41] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[12,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[12,2],round(dropRcSummer, 2 )[12,2]))
      result[vart,42] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[12,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[12,3],round(dropRcSummer, 2 )[12,3]))
      result[vart,43] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[13,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[13,1],round(dropRcSummer, 2 )[13,1]))
      result[vart,44] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[13,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[13,2],round(dropRcSummer, 2 )[13,2]))
      result[vart,45] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[13,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[13,3],round(dropRcSummer, 2 )[13,3]))
      result[vart,46] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[14,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[14,1],round(dropRcSummer, 2 )[14,1]))
      result[vart,47] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[14,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[14,2],round(dropRcSummer, 2 )[14,2]))
      result[vart,48] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[14,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[14,3],round(dropRcSummer, 2 )[14,3]))
      result[vart,49] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[15,1],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[15,1],round(dropRcSummer, 2 )[15,1]))
      result[vart,50] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[15,2],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[15,2],round(dropRcSummer, 2 )[15,2]))
      result[vart,51] <- ifelse(is.na(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]),round(ACESummer, 2 )[15,3],
                                ifelse(fitcompare$comparison[which(fitcompare$AIC==min(fitcompare$AIC))]=='AEfit',round(AESummer, 2 )[15,3],round(dropRcSummer, 2 )[15,3]))

  }, error = function(e) {print(paste("error",i,sep=" "))})
}
write.csv(result,file = '/public/home/hongxming/bmi_ace/hebi_bmi_results_2.csv')

