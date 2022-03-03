## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(saeHB.panel)
data("dataPanel")

## -----------------------------------------------------------------------------
area = max(dataPanel[,2])
period = max(dataPanel[,3])
vardir = dataPanel[,4]
result=Panel(ydi~xdi1+xdi2,area=area, period=period, vardir=vardir ,iter.mcmc = 10000,thin=5,burn.in = 1000,data=dataPanel)

## -----------------------------------------------------------------------------
result$Est

## -----------------------------------------------------------------------------
result$coefficient

## -----------------------------------------------------------------------------
result$refvar

## -----------------------------------------------------------------------------
MSE_HB=result$Est$SD^2
summary(MSE_HB)

## -----------------------------------------------------------------------------
RSE_HB=sqrt(MSE_HB)/result$Est$MEAN*100
summary(RSE_HB)

## -----------------------------------------------------------------------------
y_dir=dataPanel[,1]
y_HB=result$Est$MEAN
y=as.data.frame(cbind(y_dir,y_HB))
summary(y)
MSE_dir=dataPanel[,4]
MSE=as.data.frame(cbind(MSE_dir, MSE_HB))
summary(MSE)
RSE_dir=sqrt(MSE_dir)/y_dir*100
RSE=as.data.frame(cbind(MSE_dir, MSE_HB))
summary(RSE)

