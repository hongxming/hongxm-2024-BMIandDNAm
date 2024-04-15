rm(list = ls())
Sys.setenv(OMP_NUM_THREADS=parallel::detectCores())
library(OpenMx,lib='~/Rlib')
load('/public/home/hongxming/data_longi/deta_longi2.Rdata')
mxOption( NULL, 'Default optimizer' ,'SLSQP' )
mxOption( NULL, "mvnRelEps", value= mxOption(NULL, key="mvnRelEps")/5)
getOption('mxOptions')$mvnRelEps
CPG_all_re <- CPG_all
load('/public/home/hongxming/bmi_ace/cpg_all.Rdata')
results_bmi_longi1_acde <- read.csv('/public/home/hongxming/uni_bmi/results_bmi_longi1_o.csv')
rownames(results_bmi_longi1_acde) <- results_bmi_longi1_acde$CPG
cpg_all <- cpg_all[which(cpg_all%in%CPG_all)]
CPG_all <- cpg_all
result$A <- NA
result$Alo <- NA
result$Aup <- NA
result$CD <- NA
result$CDlo <- NA
result$CDup <- NA
result$E <- NA
result$Elo <- NA
result$Eup <- NA
lonresult <- result[,1:2]
colnames(lonresult) <- c("CPG","P_value")
for (i in 1:length(CPG_all)) {
  tryCatch({
    deta <- deta_re
    vart <- CPG_all[i]
    ###########################################################
    ###########################################################
    #set variables' parameters
    nv        <- 1                         # number of variables
    ncv       <- 2                         # number of covariates
    ntv       <- nv*2                      # number of total variables

    selVars <- c(paste0(vart,"_t1"),paste0(vart,"_t2"))
    covVars <- c("age_2_t1","age_2_t2","sex_t1","sex_t2") # You can reduce or increase the covariates as you need
    ageVars <- c("age_2_t1","age_2_t2")
    
    deta[,ageVars] = scale(deta[,ageVars])
    mzData  <- subset(deta, Zyg_t1=='MZ', c(selVars,covVars))
    dzData  <- subset(deta, Zyg_t1=='DZ', c(selVars,covVars))
    sum(is.na(mzData))
    sum(is.na(dzData))
    mzData <- na.omit(mzData)
    dzData <- na.omit(dzData)
    dim(mzData); dim(dzData)
    
    #------------ace/threshold crude estimate---------------------------
    r_MZ <- cor(mzData[,1:2],use='complete.obs')
    r_MZ
    r_DZ <- cor(dzData[,1:2],use='complete.obs')
    r_DZ
    
    #----------starting values for *mean/variance/ace* of continuous variables,sometimes maybe also suit for categorical var
    # 1 The means
    mMZ <- mean(unlist(mzData[,1:2]),na.rm=TRUE)  #only used as start values for continuous variables,ordinal variables is fixed to 0
    mDZ <- mean(unlist(dzData[,1:2]),na.rm=TRUE)  #only used as start values for continuous variables
    mMZ
    mDZ
    # 2 The covariance matrices used as start values for covariance matrices if you set "type=symm"
    cMZ <- cov(mzData[,1:2] , use='complete.obs')
    cDZ <- cov(dzData[,1:2] , use='complete.obs')
    cMZ                                           #in model write as"values=cMZ"
    cDZ
    
    # 3 The correlation matrices used as start values for covariance matrices if you set "type=stand"
    # 3.1 continuous variable
    crMZ <- r_MZ
    crDZ <- r_DZ
    crMZ                                           
    crDZ
    
    lbVas <- matrix(NA,ntv,ntv)           # lower bound for on diagonal of covariance matrix if "type=symm"
    diag(lbVas) <- 0.0001
    
    laMeMZT    <- c("m1MZT","m2MZT")          # labels for means for MZT twins
    laMeDZT    <- c("m1DZT","m2DZT")          # labels for means for DZT twins
    laVaMZT    <- c("MZ11","MZ12","MZ12","MZ22")  # labels for (co)variances for MZT twins
    laVaDZT    <- c("DZ11","DZ12","DZ12","DZ22")  # labels for (co)variances for DZT twins
    
    # 1) Specify Saturated Model with Age effects on the thresholds 
    # ---------------------------------------------------------------------------------------
    modelMZ<- mxModel("MZ",
                      # Define definition variables
                      mxMatrix( type="Full", nrow=ncv, ncol=2, free=F, label=c("data.age_2_t1","data.sex_t1",
                                                                               "data.age_2_t2","data.sex_t2"), name="MZDefVars"),
                      mxMatrix( type="Full", nrow=1, ncol=ncv, free=TRUE, values=0,labels=c("Bage","Bsex"), name="BageTH"),
                      
                      # expected mean Matrices
                      mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=mMZ, labels=laMeMZT, name="meanMZT" ), #intercepts#
                      mxAlgebra( expression= meanMZT + BageTH%*%MZDefVars, name="expMeanMZT" ), 
                      
                      mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=T, values=cMZ,labels=laVaMZT,lbound=lbVas,name="expCovMZ"),  #choose this if you want to fit submodel#
                      #mxAlgebra(expCovMZ[2,1]/sqrt(expCovMZ[1,1]%*%expCovMZ[2,2]),name="CorMZ"),
                      mxAlgebra(expression = cov2cor(expCovMZ), name="CorMZ"),
                      # Read in the data
                      mxData( observed=mzData, type="raw" ), 
                      mxExpectationNormal( covariance="expCovMZ", means="expMeanMZT", dimnames=selVars),
                      mxFitFunctionML() 
    )
    
    modelDZ<- mxModel("DZ",
                      mxMatrix( type="Full", nrow=ncv, ncol=2, free=F, label=c("data.age_2_t1","data.sex_t1",
                                                                               "data.age_2_t2","data.sex_t2"), name="DZDefVars"),
                      mxMatrix( type="Full", nrow=1, ncol=ncv, free=TRUE, values=0,labels=c("Bage","Bsex"), name="BageTH"),
                      
                      mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=mDZ, labels=laMeDZT, name="meanDZT" ), #two labels#
                      mxAlgebra( expression= meanDZT + BageTH%*%DZDefVars, name="expMeanDZT" ), 
                      mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=T, values=cDZ,labels=laVaDZT,lbound=lbVas,name="expCovDZ"), 
                      #mxAlgebra(expCovDZ[2,1]/sqrt(expCovDZ[1,1]%*%expCovDZ[2,2]),name="CorDZ"),
                      mxAlgebra(expression = cov2cor(expCovDZ), name="CorDZ"),
                      mxData( observed=dzData, type="raw"), 
                      mxExpectationNormal( covariance="expCovDZ", means="expMeanDZT", dimnames=selVars ),
                      mxFitFunctionML()
    )
    
    Conf	<- mxCI (c ('MZ.expCovMZ[2,1]','DZ.expCovDZ[2,1]',"MZ.CorMZ","DZ.CorDZ")) #,"MZ.CorMZ","DZ.CorDZ"#
    
    SatModel 	<- mxModel( "Sat", modelMZ, modelDZ, mxFitFunctionMultigroup(c('MZ.fitfunction','DZ.fitfunction')), Conf )
    
    
    SatFit	<-   mxTryHardOrdinal(SatModel, intervals=TRUE,extraTries = 15) 
    (SatSumm	<- summary(SatFit))
    
    round( SatSumm$CI[,c(2,1,3)] , 3)  #in the order: Estimate, lower bound, upper bound#
    
    
      svace <- SatFit$MZ.expCovMZ[1,1]$values/3
      laMeMZT    <- c("m1MZT","m1MZT")          # labels for means for MZT twins
      laMeDZT    <- c("m1DZT","m1DZT")          # labels for means for DZT twins
      
      modelbase <- mxModel("base",
                           mxMatrix( type="Full", nrow=1, ncol=ncv, free=TRUE, values=0,labels=c("Bage","Bsex"), name="BageTH"),
                           #mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, labels="mean", name="meanG" ),
                           # Matrices declared to store a, c, and e Path Coefficients
                           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=sqrt(svace), label="a11",lbound=0.0001, name="a" ),  #lbound=0锛燂紵#
                           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=sqrt(svace), label="d11",lbound=0.0001, name="d" ),
                           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=sqrt(svace), label="e11",lbound=0.0001, name="e" ),
                           # Matrices generated to hold A, C, and E components and total Variance
                           mxAlgebra( expression=a %*% t(a), name="A" ),
                           mxAlgebra( expression=d %*% t(d), name="D" ), 
                           mxAlgebra( expression=e %*% t(e), name="E" ),
                           mxAlgebra( expression=A+D+E, name="V" ),
                           mxAlgebra( expression=A/V, name="h2" ),
                           mxAlgebra( expression=D/V, name="d2" ),
                           mxAlgebra( expression=E/V, name="e2" ))
      
      modelMZ<- mxModel("MZ",
                        mxMatrix( type="Full", nrow=ncv, ncol=2, free=F, label=c("data.age_2_t1","data.sex_t1",
                                                                                 "data.age_2_t2","data.sex_t2"), name="DefVars"),
                        mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=mMZ, labels=laMeMZT, name="mzmeanG" ),
                        
                        mxAlgebra( expression= mzmeanG + base.BageTH%*%DefVars, name="expMean" ),
                        #mxMatrix( type="Full", nrow=nth, ncol=ntv, free=T, values=svthmz,lbound=-3, ubound=3, labels="Tmz", name="Thmz"),#thremz鈮爐hredz#
                        mxAlgebra( expression= rbind( cbind(base.A+base.D+base.E , base.A+base.D),
                                                      cbind(base.A+base.D        , base.A+base.D+base.E)), name="expCovMZ" ),
                        mxData( observed=mzData, type="raw" ), 
                        mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars ),mxFitFunctionML())
      
      modelDZ<- mxModel("DZ",
                        mxMatrix( type="Full", nrow=ncv, ncol=2, free=F, label=c("data.age_2_t1","data.sex_t1",
                                                                                 "data.age_2_t2","data.sex_t2"), name="DefVars"),
                        mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=mDZ, labels=laMeDZT, name="dzmeanG" ),
                        
                        mxAlgebra( expression= dzmeanG + base.BageTH%*%DefVars, name="expMean" ),
                        #mxMatrix( type="Full", nrow=nth, ncol=ntv, free=T, values=svthdz,lbound=-3, ubound=3, labels="Tdz", name="Thdz"),#thremz鈮爐hredz#
                        mxAlgebra( expression= rbind( cbind(base.A+base.D+base.E , 0.5%x%base.A+0.25%x%base.D),
                                                      cbind(0.5%x%base.A+0.25%x%base.D  , base.A+base.D+base.E)), name="expCovDZ" ),
                        mxData( observed=dzData, type="raw" ),
                        mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars),mxFitFunctionML())
      
      # Create Confidence Interval Objects
      
      Conf		<- mxCI (c("base.A[1,1]","base.D[1,1]","base.E[1,1]","base.h2[1,1]","base.d2[1,1]","base.e2[1,1]" ))
      
      AdeModel 	<- mxModel( "ADE",modelbase, modelMZ, modelDZ, 
                            mxFitFunctionMultigroup(c('MZ.fitfunction','DZ.fitfunction')),Conf)
      
      #--------------------------------------------------------------------------
      # 4) RUN AceModel
      #AceFit   <- mxRun(AceModel, intervals=TRUE)
      AdeFit   <- mxTryHardOrdinal(AdeModel, intervals=TRUE,extraTries = 15)
      (AdeSumm  <- summary(AdeFit,verbose=T))
      # aderesult <- list(AdeSumm$CIdetail,AdeSumm$CI,subresult)
      AdeSumm$CI[4,1] <- AdeSumm$CIdetail[7,3]
      AdeSumm$CI[5,1] <- AdeSumm$CIdetail[9,3]
      AdeSumm$CI[6,1] <- AdeSumm$CIdetail[11,3]
      AdeSumm$CI[4,3] <- AdeSumm$CIdetail[8,3]
      AdeSumm$CI[5,3] <- AdeSumm$CIdetail[10,3]
      AdeSumm$CI[6,3] <- AdeSumm$CIdetail[12,3]
      
      lon <- omxSetParameters(model=AdeFit,labels="e11",free=F,values=sqrt(results_bmi_longi1_acde[vart,'e2']),name="lon")
      lon <- omxSetParameters(model=lon,labels="a11",free=F,values=sqrt(results_bmi_longi1_acde[vart,'a2']),name="lon")
      lonfit <- mxTryHardOrdinal(lon,intervals=T,extraTries = 15)
      loncompare <- mxCompare(base=AdeFit,comparison=lonfit)
      mxCompare (SatFit,AdeFit)
    lonresult[vart,'P_value'] <- loncompare[2,'p']
  }, error = function(e) {print(paste("error",i,sep=" "))})
}
write.csv(lonresult,file = '/public/home/hongxming/bmi_ace/lonresults_bmi_acde.csv')

