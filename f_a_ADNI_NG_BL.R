#!/usr/bin/env Rscript
################################################################################
## ADNI-NG-BL workflow - PART A::MULTIVARIATE LINEAR MODELS ####################
############################################################ MAFF V001 #########
library(Biobase)
library(limma)
################################################################################
# input:
#      + NG metabolite matrix (Excel workbook file)
#      + ADNI main covariates file (csv file)
# out:
#      + linear model with several AD contrasts (Excel workbook file)
################################################################################
f_input2eset <- function(dat){
  metFile<-"mt_ready_ADNI1_GO2_Baseline_Nightingale.xlsx"
  mx<-readxl::read_excel(path=paste0(getwd(),"/", metFile), sheet=1, na = "")
  meta<-readxl::read_excel(path=paste0(getwd(), "/", metFile), sheet=2, na = "")
  mx<-cbind(RID=meta[, 27], mx[, 2:ncol(mx)])
  meta<-cbind(meta[, 1], meta[, 26:27], meta[, 2:25], meta[, 28:70], meta[, 72])
  meta <-meta[order(meta$RID),]
  mx<-mx[order(mx$RID),]
  rownames(mx)<-mx[, 1]
  #
  metAnno<-readxl::read_excel(path=paste0(getwd(), "/", metFile),sheet=3, na="")
  adniMerge<-read.csv(file=paste0(getwd(),"/ADNIMERGE.csv"), 
                      header=TRUE, na.strings=c("NA"," ", "", ""))
  adniMerge<-adniMerge[which(adniMerge$RID %in% mx$RID),]
  adniMerge<-subset(adniMerge, VISCODE == "bl")
  adniMerge<-adniMerge[order(adniMerge$RID),]
  #
  met<- t(mx)
  colnames(met)<-met[1, ]
  met<-met[c(2:nrow(met)), ]
  colsData<-cbind(adniMerge[, c(1,4,9:15)], adniMerge[, grep("_bl$", 
                                          colnames(adniMerge))], meta[, -c(3)])
  row.names(colsData)<-colsData[, 1]
  colsData<-colsData[, -1]
  print(paste0("colnames(met) == row.names(colsData): ", 
                                    all(colnames(met) == row.names(colsData))))
  phenoData<-new("AnnotatedDataFrame", data=colsData)
  eset<-Biobase::ExpressionSet(assayData=as.matrix(met), phenoData=phenoData)
  return(eset)
}

meset<-f_input2eset(dat)
print(meset)
#------------------------------------------------------------------------------#

print(table(meset$DX_bl))
meset$DX_bl<-factor(meset$DX_bl)
meset$PTGENDER<-factor(meset$PTGENDER)
meset$PTRACCAT<-factor(meset$PTRACCAT)
meset$PTRACCAT<-gsub("More than one","Other",meset$PTRACCAT)
meset$PTRACCAT<-gsub("Am Indian/Alaskan","Other",meset$PTRACCAT)
meset$PTRACCAT<-gsub("Hawaiian/Other PI","Other",meset$PTRACCAT)
meset$PTRACCAT<-gsub("Unknown","Other",meset$PTRACCAT)

meset$PTEDUCAT<-ifelse(meset$PTEDUCAT <5, "0-5", ifelse(meset$PTEDUCAT <9,"5-9", 
                                  ifelse(meset$PTEDUCAT <14, "10-14", "15-20")))
meset$PTEDUCAT<-as.factor(meset$PTEDUCAT)

design <- model.matrix(~0 + DX_bl + AGE + PTGENDER + PTRACCAT + 
                                              PTEDUCAT + Med.Statin, data=meset)
colnames(design) = gsub("DX_bl","",colnames(design))
colnames(design) = make.names(colnames(design))
x <- c("AD-CN","AD-SMC","AD-EMCI","AD-LMCI")
contrast = makeContrasts(contrasts=x,levels=design)
fit1 <- lmFit(exprs(meset), design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)

# make out directory
check_create_dir <- function(the_dir) {
  if (!dir.exists(the_dir)) {
    dir.create(the_dir, recursive = TRUE) }
}
the_dir <- "./results"
check_create_dir(the_dir)

l.res.ad.cn = topTable(fit3, coef=c("AD-CN"), n=Inf)
l.res.ad.smc = topTable(fit3, coef=c("AD-SMC"), n=Inf)
l.res.ad.emci = topTable(fit3, coef=c("AD-EMCI"), n=Inf)
l.res.ad.lmci = topTable(fit3, coef=c("AD-LMCI"), n=Inf)


## add Cohen's effect size (d)
# fit$coefficients / sqrt(fit$s2.post)
# logFC / sqrt(fit3$s2.post) FOR each contrast

s2<- data.frame(MET=names(fit3$s2.post), P.STDEV=fit3$s2.post)
s2<- s2[sort(s2$MET),]

# ADvsCN
ad.cn<-cbind(l.res.ad.cn[sort(row.names(l.res.ad.cn)),], s2)
ad.cn$CohenD<-ad.cn$logFC / sqrt(ad.cn$P.STDEV)
ad.cn<-ad.cn[, -c(7,8)]

# ADvsSMC
ad.smc<-cbind(l.res.ad.smc[sort(row.names(l.res.ad.smc)),], s2)
ad.smc$CohenD<-ad.smc$logFC / sqrt(ad.smc$P.STDEV)
ad.smc<-ad.smc[, -c(7,8)]

# ADvsEMCI
ad.emci<-cbind(l.res.ad.emci[sort(row.names(l.res.ad.emci)),], s2)
ad.emci$CohenD<-ad.emci$logFC / sqrt(ad.emci$P.STDEV)
ad.emci<-ad.emci[, -c(7,8)]

# ADvsLMCI
ad.lmci<-cbind(l.res.ad.lmci[sort(row.names(l.res.ad.lmci)),], s2)
ad.lmci$CohenD<-ad.lmci$logFC / sqrt(ad.lmci$P.STDEV)
ad.lmci<-ad.lmci[, -c(7,8)]

# write to Excel
outdat = list(ADvsCN=ad.cn, ADvsSMC=ad.smc, 
                ADvsEMCI=ad.emci, ADvsLMCI=ad.lmci)
openxlsx::write.xlsx(outdat, file='./results/ADNI-NG_model_1.xlsx', 
                                                                 row.names=TRUE)
#------------------------------------------------------------------------------#

# TODO -->> to mirror Jun's UKBB NG study
# FROM ALEJO notes:

# Model 1 = Prevalent all cause dementia ~ meta + age + sex + BMI + 
# batchnightingale + spectrometer number + fasting time + race. 
# We will use the uncorrected dataset from ADNI, and the regression 
# will be logistic.

# Model 2 = model 1 + smoking + alcohol + education. We will use the 
# medication-adjusted dataset from ADNI. This will take care of some further 
# covariants that Jun used in UK Biobank to control medication usage. 
# The covariants that Jun used in UK Biobank were "metabolic medications 
# (Jun to send) + CNS medications (4 digits ATC, n = 18) + 
# depression medications + proton pump inhibitor". Jun will check whether 
# the scale of the resulting beta values from ADNI (when correcting medication 
# in two steps) are going to be in the same scale as the beta values from UKBB

# Model 3 = model 2 + APOEe4. Once model 2 has been created, extending this 
# model into model 3 should not be difficult
#------------------------------------------------------------------------------#