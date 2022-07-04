# hdim.ui
hdim.ui - R package for sensitivity analysis in high-dimensional causal inference and missing outcome data contexts

Variables:
(1) age, (2) education (less than high school, high
                        school, greater than high school), (3) white (yes/no), (4) hispanic (yes/no), (5) foreign (yes/no),
(6) alcohol use, (7) married (yes/no), (8) liver birth order (one, two, greater than 2), (9) number
of prenatal visits.3


```
library(foreign)
data=read.dta("C_2010_JOE-dataRS5K.dta")
length(names(data))

categoricalCovariates=c("dmar","mwhite","mhispan","foreignb","alcohol","dmeduc","dlivord")
#quantitative
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

withInteractions= model.matrix(~ .^30, data = newData[,c(categoricalCovariates,continuousCovariates)],raw = TRUE)[,-1]
higherorder=poly(cbind(newData[,"dmage"], newData[,"nprevist"]),degree=20,raw = TRUE)


length(names(newData))

# remove dmage and nprevist interaction bewteen them since it is included in withInteractions
higherorder=higherorder[,-c(1,4,5)]

covariates=cbind(withInteractions,higherorder)
ncol(covariates)
Y=newData[,"dbirwt"]
T=ifelse(newData[,"T"]>=1,1,0)

ui_real=ui.causal(X=covariates,Y=(Y),T=(T),rho0=c(-0.15,0.15),rho1=c(-0.15,0.15),
                  # subset="double",
                  subset = "refit", sigma_correction = 'new',
                   param="ACE", rho.plotrange = c(-0.25, 0.25),regularization_alpha = 1)


plot(ui_real,ylab="ACE")
ui_y0t1=ui.causal(X=covariates,Y=(Y),T=(T),rho0=c(-0.3,0.1),rho.plotrange =c(-0.3,0.1),gridn=41,
                  # subset="double",
                  subset = "refit", sigma_correction = 'new',
                  param="Y0T1",regularization_alpha = 1)
ui_y1t0=ui.causal(X=covariates,Y=(Y),T=(T),rho1=c(-0.3,0.1),rho.plotrange =c(-0.3,0.1),gridn=41,
                  # subset="double",
                  subset = "refit", sigma_correction = 'new',
                  param="Y1T0",regularization_alpha = 1)



# 3130
Y1T1_mean=mean(Y[T==1])
# 3410
Y0T0_mean=mean(Y[T==0])
plot(ui_y0t1,infer = FALSE, lower=Y1T1_mean,
     upper=Y0T0_mean,xlab=expression(rho[0]),ylab=expression(E(Y(0)*"|"*T*"="*1)))
plot(ui_y1t0,infer = FALSE, lower=Y1T1_mean,upper=Y0T0_mean
     ,xlab=expression(rho[1]),ylab=expression(E(Y(1)*"|"*T*"="*0)))

#rh0 between -0.25 and 0.05
ui_y0t1$DR$coef
#rh1 between -0.2 and 0.06
ui_y1t0$DR$coef

# ui_real2=ui(X=covariates,Y=(Y),T=(T),rho0=c(-0.25,0.05),rho1=c(-0.2,0.05),subset="refit",param="ACE",regularization_alpha = 1)
ui_real2=ui.causal(X=covariates,Y=(Y),T=(T),rho0=c(-0.25,0.05),rho1=c(-0.2,0.06),
                   # subset="double",
                   subset = "refit", sigma_correction = 'new',
                   param="ACE",regularization_alpha = 1)
plot(ui_real2)

# (-332.513, 50.194)
print(ui_real2)
# -279.2004
min(ui_real2$DR$coef)
# -3.118664
max(ui_real2$DR$coef)


# -219.5713
ui_real2$DR$coef["0","0"]
# -272.8838
ui_real2$DR$ci["0","0",1]
# -166.2589
ui_real2$DR$ci["0","0",2]
```
