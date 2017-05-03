########################################################################################################
# Panel data analysis of ICMS-E effects in Brazil # Sept 2016 # nils.droste@ufz.de
# research question: is there an effect of ICMS-E on nature conservation area (%) among states in Brazil?
########################################################################################################

#### set a working directory ####

library(plm)
library(lmtest)
library(sandwich)
library(stargazer)
library(punitroots)

# setwd("/home/HereGoesYourWD/")
df <- read.csv(paste(getwd(), "/paneldataEFT-BR.csv", sep=""), sep=",", dec=".")
df <- pdata.frame(df, index = c("ID", "year"))

# Info on variables --------------------------------------------------------

#Ys
#tot<-total protected area share of state territory in %
#fed<-federal protected area share of state territory in %
#sta<-state protected area share of state territory in %
#mun<-municipal protected area share of state territory in %
#tPI<-total protected area (protecao integral) share of state territory in %
#fPI<-federal protected area (protecao integral) share of state territory in %
#sPI<-state protected area (protecao integral) share of state territory in %
#mPI<-municipal protected area (protecao integral) share of state territory in %
#tUS<-total protected area (uso sustentavel) share of state territory in %
#fUS<-federal protected area (uso sustentavel) share of state territory in %
#sUS<-state protected area (uso sustentavel) share of state territory in %
#mUS<-municipal protected area (uso sustentavel) share of state territory in %

#Xs
#ID<-states of brazil
#year<-year
#icms_e<- 1 if ecological fiscal transfers enactment in state i year t, 0 if not
#icms_e1<- 1 if ecological fiscal transfers first legislative decision in state i year t, 0 if not
#agr<-share of valued added by agriculture of valued added by economic activity (agr+ind+ser) in % (constant prices R$2000) in thousands
#ind<-share of valued added by industry of valued added by economic activity (agr+ind+ser) in % (constant prices R$2000)in thousands
#ser<-share of valued added by service of valued added by economic activity (agr+ind+ser) in % (constant prices R$2000)in thousands
#pop<-population density cap/km²
#inc<-GDP per capita in constant prices R$2000 in thousands
#arpa<-program for protected areas in the amazon 1 if in force in state i year t, 0 if not
#ama<- amazon biome 1 if major share of state i territory
#cer<- cerrado biome 1 if major share of state i territory
#caa<- caatinga biome 1 if major share of state i territory
#mat<- mata atlantica biome 1 if major share of state i territory
#pan<- pantanal biome 1 if major share of state i territory
#pam<- pampa biome 1 if major share of state i territory


# start -------------------------------------------------------------------

#check panel df and write an html summary table
str(df)
summary(df)
stargazer(df[c("tot","fed","sta","mun","icms_e","agr","ind","ser","pop","inc")], type = "html", title="Descriptive statistics", digits=1, out="summary2.html",covariate.labels=c("total protected area share of state territory in per cent (tot)","federal protected area share of state territory in per cent (fed)","state protected area share of state territory in per cent (sta)","municipal protected area share of state territory in per cent (mun)","ICMS-E dummy (icmse)","share of valued added by agriculture in per cent (agr)","share of valued added by industry in per cent (ser)","share of valued added by service in per cent (ind)","population density cap/km² (pop)","GDP per capita, R$ in thousands (inc)"))

#constant = 0.5 times the minimal observed value #required for log transformations of variables /w 0
c<-min(df$mun[which(df$mun > 0)])*0.5
c2<-min(df$sta[which(df$sta > 0)])*0.5

#testing for unit root vs stationarity with Im-Pesaran-Shin test
ptestdat<-pdata.frame(as.data.frame(cbind(ID=df$ID,year=df$year,ltot=log(df$tot),lfed=log(df$fed),lstat=log(df$sta+c2),lmun=log(df$mun+c),lagr=log(df$agr),lind=log(df$ind),lpop=log(df$pop),linc=log(df$inc))), index=c("ID","year"))[,3:10]
pt<-purtest(ptestdat, test = "ips", exo="trend", Hcons=T, lags="AIC", pmax = 10); pt # rejects unit root # but this assumes cross-sectional independence
summary(pt)

# determining order of integration through a cross-sectionally augmented IPS Tests for single time series
# basically all test on time series, analysed individually, do not reject the Null of a Unit-Root until fifferenced
cipstest(diff(ptestdat$ltot), lags = pt$idres$ltot$lags, type = "trend", model = "dmg") # OoI (order of integration): 1
cipstest(diff(ptestdat$lfed), lags = pt$idres$lfed$lags, type = "trend", model = "dmg") # OoI (order of integration): 1
cipstest(diff(diff(diff(ptestdat$lstat))), lags = pt$idres$lstat$lags, type = "trend", model = "dmg") # not determinable
cipstest(diff(diff(ptestdat$lmun)), lags = pt$idres$lmun$lags, type = "trend", model = "dmg") # OoI (order of integration): 2
cipstest(diff(ptestdat$lagr), lags = pt$idres$lagr$lags, type = "trend", model = "dmg") # OoI (order of integration): 1
cipstest(diff(ptestdat$lind), lags = pt$idres$lind$lags, type = "trend", model = "dmg") # OoI (order of integration): 1
cipstest(diff(ptestdat$lpop), lags = pt$idres$lpop$lags, type = "trend", model = "dmg") # OoI (order of integration): 1
cipstest(diff(ptestdat$linc), lags = pt$idres$linc$lags, type = "trend", model = "dmg") # OoI (order of integration): 1

# the CADF test on the panel without differencing on the other hand rejects the Null quite clearly
Demetrescu <- pCADFtest(ptestdat, covariates = NULL, crosscorr=.05, type="trend", max.lag.y=10, criterion="AIC"); Demetrescu
summary(Demetrescu)

# all PA regressions ------------------------------------------------------

#oneway random effects (model 1)
m1.tot<- plm(log(tot)~icms_e+log(agr)+log(ind)+log(pop)+log(inc)+arpa, data=df, index=c("ID", "year"), model="random", effect="individual", random.method = "swar")
summary(m1.tot)
pbgtest(m1.tot) # null is that there is no serial correlation, p < 0.5 <- serial correlation
pcdtest(m1.tot) # null is that there is no cross-sectional dependence , p < 0.5 <- cross-sectional dependence
m1.tot.rob<-coeftest(m1.tot, vcov=function(x) vcovSCC(x, type="HC3", maxlag=10)); m1.tot.rob

# m1 + biomes (model 2)
m2.tot <- plm(log(tot)~icms_e+log(agr)+log(ind)+log(pop)+log(inc)+arpa+ama+cer+caa+mat+pan+pam,df, index=c("ID", "year"), model="random", effect="individual", random.method = "swar")
summary(m2.tot)
pbgtest(m2.tot) # null is that there is no serial correlation, p < 0.5 <- serial correlation
pcdtest(m2.tot) # null is that there is no cross-sectional dependence , p < 0.5 <- cross-sectional dependence
m2.tot.rob<-coeftest(m2.tot, .vcov=function(x) vcovSCC(x, type="HC3", maxlag=2)); m2.tot.rob

# m2 + detrend (model 3)
m3.tot<-plm(log(tot)~icms_e+log(agr)+log(ind)+log(pop)+log(inc)+arpa+ama+cer+caa+mat+pan+pam+icms_e*as.numeric(year), data=df, index=c("ID", "year"), model="random", effect="individual")
summary(m3.tot)
pbgtest(m3.tot) # null is that there is no serial correlation, p < 0.5 <- serial correlation
pcdtest(m3.tot) # null is that there is no cross-sectional dependence , p < 0.5 <- cross-sectional dependence
m3.tot.rob<-coeftest(m3.tot, vcov=function(x) vcovSCC(x, type="HC3", maxlag=2)); m3.tot.rob

# m3 + interaction (model 4)
m4.tot<-plm(log(tot)~icms_e+log(agr)+log(ind)+log(pop)+log(inc)+arpa+ama+cer+caa+mat+pan+pam+icms_e*as.numeric(year)+icms_e*log(agr)+icms_e*log(ind)+icms_e*log(pop)+icms_e*log(inc), data=df, index=c("ID", "year"), model="random", effect="individual")
summary(m4.tot)
pbgtest(m4.tot) # null is that there is no serial correlation, p < 0.5 <- serial correlation
pcdtest(m4.tot) # null is that there is no cross-sectional dependence , p < 0.5 <- cross-sectional dependence
m4.tot.rob<-coeftest(m4.tot, vcov=function(x) vcovSCC(x, type="HC3", maxlag=2)); m4.tot.rob

#analytical plots
par(mfrow=c(2,2))
plot(density(resid(m1.tot)),main = "residual density")
fitted.m1.tot=(m1.tot$model[[1]] - m1.tot$residuals)
plot(resid(m1.tot)~(fitted.m1.tot), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.m1.tot, residuals(m1.tot)), col="red")
qqnorm(m1.tot$resid)
qqline(m1.tot$resid, col="red")
lev = hat(model.matrix(m1.tot))
plot(lev, main="leverage")
#df[lev>0.4,]

plot(density(resid(m2.tot)),main = "residual density")
fitted.m2.tot=(m2.tot$model[[1]] - m2.tot$residuals)
plot(resid(m2.tot)~(fitted.m2.tot), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.m2.tot, residuals(m2.tot)), col="red")
qqnorm(m2.tot$resid)
qqline(m2.tot$resid, col="red")
lev = hat(model.matrix(m2.tot))
plot(lev, main="leverage")
#df[lev>0.4,]

plot(density(resid(m3.tot)),main = "residual density")
fitted.m3.tot=(m3.tot$model[[1]] - m3.tot$residuals)
plot(resid(m1.tot)~(fitted.m3.tot), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.m3.tot, residuals(m3.tot)), col="red")
qqnorm(m3.tot$resid)
qqline(m3.tot$resid, col="red")
lev = hat(model.matrix(m3.tot))
plot(lev, main="leverage")
#df[lev>0.4,]

plot(density(resid(m4.tot)),main = "residual density")
fitted.m4.tot=(m4.tot$model[[1]] - m4.tot$residuals)
plot(resid(m4.tot)~(fitted.m4.tot), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.m4.tot, residuals(m4.tot)), col="red")
qqnorm(m4.tot$resid)
qqline(m4.tot$resid, col="red")
lev = hat(model.matrix(m4.tot))
plot(lev, main="leverage")
#df[lev>0.4,]
#dev.off()

#output table
rob.se<-list(m1.tot.rob[,2],m2.tot.rob[,2],m3.tot.rob[,2],m4.tot.rob[,2])
rob.p<-list(m1.tot.rob[,4],m2.tot.rob[,4],m3.tot.rob[,4],m4.tot.rob[,4])
stargazer(m1.tot, m2.tot, m3.tot, m4.tot, dep.var.labels="ln of protected area share in percent of total area", p=rob.p, se=rob.se, type="html", out="table_fed.htm")

# municipal regressions ----------------------------------------------------
library(plm)
library(lmtest)
library(sandwich)
library(stargazer)

#basic random effects regression (model 1)
m1.mun <- plm(log(mun+c)~icms_e+log(agr)+log(ind)+log(pop)+log(inc)+log(fed)+log(sta+c2)+arpa, data=df, index=c("ID", "year"), model="random", effect="individual")
summary(m1.mun)
pbgtest(m1.mun) # null is that there is no serial correlation, p < 0.5 <- serial correlation
pcdtest(m1.mun) # null is that there is no cross-sectional dependence , p < 0.5 <- cross-sectional dependence
m1.mun.rob<-coeftest(m1.mun, vcov=function(x) vcovSCC(x, type="HC3", maxlag=2)); m1.mun.rob

#random effects regression with biomes (model 2)
m2.mun <- plm(log(mun+c)~icms_e+log(agr)+log(ind)+log(pop)+log(inc)+log(fed)+log(sta+c2)+arpa+ama+cer+caa+mat+pan+pam,df, index=c("ID", "year"), model="random", effect="ind")
summary(m2.mun)
pbgtest(m2.mun) # null is that there is no serial correlation, p < 0.5 <- serial correlation
pcdtest(m2.mun) # null is that there is no cross-sectional dependence , p < 0.5 <- cross-sectional dependence
m2.mun.rob<-coeftest(m2.mun, vcov=function(x) vcovSCC(x, type="HC3", maxlag=2)); m2.mun.rob

#oneway+detrend
m3.mun<-plm(log(mun+c)~icms_e+log(agr)+log(ind)+log(pop)+log(inc)+log(fed)+log(sta+c2)+arpa+ama+cer+caa+mat+pan+pam+icms_e*as.numeric(year), data=df, index=c("ID", "year"), model="random", effect="individual")
summary(m3.mun)
pbgtest(m3.mun) # null is that there is no serial correlation, p < 0.5 <- serial correlation
pcdtest(m3.mun) # null is that there is no cross-sectional dependence , p < 0.5 <- cross-sectional dependence
m3.mun.rob<-coeftest(m3.mun, vcov=function(x) vcovSCC(x, type="HC3", maxlag=2)); m3.mun.rob

#oneway+detrend+interaction
m4.mun<-plm(log(mun+c)~icms_e+log(agr)+log(ind)+log(pop)+log(inc)+log(fed)+log(sta+c2)+arpa+ama+cer+caa+mat+pan+pam+icms_e*as.numeric(year)+icms_e:log(agr)+icms_e:log(ind)+icms_e:log(pop)+icms_e:log(inc)+icms_e:log(fed)+icms_e:log(sta+c2), data=df, index=c("ID", "year"), model="random", effect="individual")
summary(m4.mun)
pbgtest(m4.mun) # null is that there is no serial correlation, p < 0.5 <- serial correlation
pcdtest(m4.mun) # null is that there is no cross-sectional dependence , p < 0.5 <- cross-sectional dependence
m4.mun.rob<-coeftest(m4.mun, vcov=function(x) vcovSCC(x, type="HC3", maxlag=2)); m4.mun.rob

#analytical plots
par(mfrow=c(2,2))
plot(density(resid(m1.mun)),main = "residual density")
fitted.m1.mun=(m1.mun$model[[1]] - m1.mun$residuals)
plot(resid(m1.mun)~(fitted.m1.mun), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.m1.mun, residuals(m1.mun)), col="red")
qqnorm(m1.mun$resid)
qqline(m1.mun$resid, col="red")
lev = hat(model.matrix(m1.mun))
plot(lev, main="leverage")

plot(density(resid(m2.mun)),main = "residual density")
fitted.m2.mun=(m2.mun$model[[1]] - m2.mun$residuals)
plot(resid(m2.mun)~(fitted.m2.mun), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.m2.mun, residuals(m2.mun)), col="red")
qqnorm(m2.mun$resid)
qqline(m2.mun$resid, col="red")
lev = hat(model.matrix(m2.mun))
plot(lev, main="leverage")
#df[lev>0.4,]

plot(density(resid(m3.mun)),main = "residual density")
fitted.m3.mun=(m3.mun$model[[1]] - m3.mun$residuals)
plot(resid(m3.mun)~(fitted.m3.mun), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.m3.mun, residuals(m3.mun)), col="red")
qqnorm(m3.mun$resid)
qqline(m3.mun$resid, col="red")
lev = hat(model.matrix(m3.mun)); 
plot(lev, main="leverage")
#df[lev>0.3,]

plot(density(resid(m4.mun)),main = "residual density")
fitted.m4.mun=(m4.mun$model[[1]] - m4.mun$residuals)
plot(resid(m4.mun)~(fitted.m4.mun), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.m4.mun, residuals(m4.mun)), col="red")
qqnorm(m4.mun$resid)
qqline(m4.mun$resid, col="red")
lev = hat(model.matrix(m4.mun))
plot(lev, main="leverage")
#df[lev>0.3,]

#output table
rob.se<-list(m1.mun.rob[,2],m2.mun.rob[,2],m3.mun.rob[,2],m4.mun.rob[,2])
rob.p<-list(m1.mun.rob[,4],m2.mun.rob[,4],m3.mun.rob[,4],m4.mun.rob[,4])
stargazer(m1.mun, m2.mun, m3.mun, m4.mun,dep.var.labels="ln of municipal protected area share in percent of total area", p=rob.p, se=rob.se, type="html", out="table_mun.html")

# one could compare consistency of random effects models via Hausman test
# phtest(fixed.mod, random.mod)