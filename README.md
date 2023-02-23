# hdim.ui
hdim.ui - R package for sensitivity analysis in high-dimensional causal inference and missing outcome data contexts

Built using code from the R-package **ui** available on CRAN.

The code below replicates the case study in the manuscript "Valid causal inference with unobserved confounding in high-dimensional settings" by Niloofar Moosavi, Tetiana Gorbach and Xavier de Luna (2023).

Please note that the link to the data provided in the code is currently pointing to a specific location. In order to access the data, you will need to download it from this link: https://github.com/mdcattaneo/replication-C_2010_JOE.

Variables:
(1) age, (2) education (less than high school, high
                        school, greater than high school), (3) white (yes/no), (4) hispanic (yes/no), (5) foreign (yes/no),
(6) alcohol use, (7) married (yes/no), (8) liver birth order (one, two, greater than 2), (9) number
of prenatal visits.3


```

install.packages(c('devtools','hdm','glmnet','foreign'))

library(devtools)
install_github('https://github.com/stat4reg/hdim.ui.git')


library(hdm)
library(glmnet)
library(hdim.ui)

# prepare data####
library(foreign)
data=read.dta("C_2010_JOE-dataRS5K.dta")
length(names(data))

# List of variables
# (1) age, (2) education (less than high school, high
#                         school, greater than high school), (3) white (yes/no), (4) hispanic (yes/no), (5) foreign (yes/no),
# (6) alcohol use, (7) married (yes/no), (8) liver birth order (one, two, greater than 2), (9) number
# of prenatal visits.3



categoricalCovariates=c("dmar","mwhite","mhispan","foreignb","alcohol","dmeduc","dlivord")
continuousCovariates=c("dmage","nprevist")

newData= data[,c("T","dbirwt",categoricalCovariates,continuousCovariates)]

newData$dmeduc[data$dmeduc <12] = "lthighschool"
newData$dmeduc[data$dmeduc == 12] = "highschool"
newData$dmeduc[data$dmeduc >12] = "mthighschool"
newData$dmeduc=as.factor(newData$dmeduc)

newData$dlivord[data$dlivord ==1] = "one"
newData$dlivord[data$dlivord ==2] = "two"
newData$dlivord[data$dlivord >2] = "gttwo"
newData$dlivord=as.factor(newData$dlivord)

length(names(newData))
withInteractions= model.matrix(~ .^2, data = newData[,c(categoricalCovariates,continuousCovariates)])[,-1]
higherorder=poly(cbind(newData[,"dmage"], newData[,"nprevist"]),degree=3)


length(names(newData))

# remove dmage and nprevist interaction bewteen them since it is included in withInteractions
higherorder=higherorder[,-c(1,4,5)]

covariates=cbind(withInteractions,higherorder)
ncol(covariates)
Y=newData[,"dbirwt"]
T=ifelse(newData[,"T"]>=1,1,0)


##########################
ui_y0t1=ui.causal(X=covariates,Y=(Y),T=(T),rho0=c(-0.3,0.1),rho.plotrange =c(-0.4,0.4),
                  # subset="double",
                  subset = "refit", sigma_correction = 'new',
                  param="Y0T1",regularization_alpha = 1)
ui_y1t0=ui.causal(X=covariates,Y=(Y),T=(T),rho1=c(-0.3,0.1),rho.plotrange =c(-0.4,0.4),
                  # subset="double",
                  subset = "refit", sigma_correction = 'new',
                  param="Y1T0",regularization_alpha = 1)



# 3132.952
Y1T1_mean=mean(Y[T==1])
# 3411.788
Y0T0_mean=mean(Y[T==0])



par(mfrow=c(1,2))

plot.ui(ui_y0t1,infer = FALSE, lower=Y1T1_mean,
     upper=Y0T0_mean,xlab=expression(rho[0]),ylab=expression(E(Y(0)*"|"*T*"="*1)))
plot.ui(ui_y1t0,infer = FALSE, lower=Y1T1_mean,upper=Y0T0_mean
     ,xlab=expression(rho[1]),ylab=expression(E(Y(1)*"|"*T*"="*0)))

############################
#rh0 between -0.25 and 0.05
ui_y0t1$DR$coef
#rh1 between -0.2 and 0.06
ui_y1t0$DR$coef

ui_y0=ui.causal(X=covariates,Y=(Y),T=(T),rho0=c(-0.25,0.05),rho.plotrange =c(-0.5,0.5),
                  # subset="double",
                  subset = "refit", sigma_correction = 'new',
                  param="Y0",regularization_alpha = 1)
ui_y1=ui.causal(X=covariates,Y=(Y),T=(T),rho1=c(-0.2,0.06),rho.plotrange =c(-0.5,0.5),
                  # subset="double",
                  subset = "refit", sigma_correction = 'new',
                  param="Y1",regularization_alpha = 1)



plot.ui(ui_y0,xlab=expression(rho[0]),ylab=expression(E(Y(0))))
plot.ui(ui_y1,xlab=expression(rho[1]),ylab=expression(E(Y(1))))


```
