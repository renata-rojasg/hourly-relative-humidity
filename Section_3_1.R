rm(list = ls())
library(forecast)
library(tidyverse)
library(BTSR)
source("ubxiiarma.fit.r")
source("best.ubxiiarma.r")
######################
## Data preparation ##
######################
data <- readr::read_delim("combined_hourly_data.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(timestamp=as.POSIXct(timestamp, tz="GMT",
                              origin="1970-01-01 00:00:00"),
         hum=hum/100)

data<-data[830:1870,] # only winter
n<-round(dim(data)[1]*.8)

#########################
## Train and test sets ##
#########################
datatrain<-cbind(data[1:n,])  
colnames(datatrain)<-paste0(colnames(datatrain),"_train")

datatest<-cbind(data[(n+1):(dim(data)[1]),])  
colnames(datatest)<-paste0(colnames(datatest),"_test")

datatrain$RSSI_03_train[is.na(datatrain$RSSI_03_train)]<-
  mean(na.omit(datatrain$RSSI_03_train))
datatrain$RSSI_04_train[is.na(datatrain$RSSI_04_train)]<-
  mean(na.omit(datatrain$RSSI_04_train))
datatrain$RSSI_05_train[is.na(datatrain$RSSI_05_train)]<-
  mean(na.omit(datatrain$RSSI_05_train))
datatrain$RSSI_06_train[is.na(datatrain$RSSI_06_train)]<-
  mean(na.omit(datatrain$RSSI_06_train))
datatrain$RSSI_08_train[is.na(datatrain$RSSI_08_train)]<-
  mean(na.omit(datatrain$RSSI_08_train))

suppressMessages(attach(datatrain))
suppressMessages(attach(datatest))
suppressMessages(attach(data))
X<-as.matrix(datatrain[,2])
Xtest<-as.matrix(datatest[,2])
X0<-rbind(X,Xtest)
nX<-dim(X)[2]
######################
## stacionary tests ##
######################
truncation<-round(12*(n/100)^(.25)) # Schwert rule
adf.level<-tseries::adf.test(hum_train,k=truncation) # stacionary
adf.diff<-tseries::adf.test(diff(hum_train),k=truncation) # stacionary
kpss.level<-tseries::kpss.test(hum_train,lshort = F) # non-stacionary
kpss.diff<-tseries::kpss.test(diff(hum_train),lshort = F) # stacionary

table.stationarity<-data.frame(
  Series=c("In Level","1st difference"),
  ADF=c(adf.level$statistic,adf.diff$statistic),
  `p-value ADF`=c(adf.level$p.value,adf.diff$p.value),
  KPSS=c(kpss.level$statistic,kpss.diff$statistic),
  `p-value KPSS`=c(kpss.level$p.value,kpss.diff$p.value)
)

########################
## Fitting the models ##
########################
# fiting the ARIMA
a01<-auto.arima(hum_train, allowdrift = F)
new1<-Arima(hum_test,model=a01) #one-step-ahead
a02<-auto.arima(hum_train, xreg = X,allowdrift = F)
new2<-Arima(hum_test,xreg = Xtest,model=a02) #one-step-ahead

# fiting the unit ARMA
quant<-.5
order<-matrix(NA,16,8)
cont<-1
for(i in 0:3){
  for(j in 0:3){
    barma<-summary(BARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,
                               report=F))
    karma1<-suppressWarnings(KARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,
                                         rho=quant,
                                         control = list(method="Nelder-Mead",stopcr=1e-2),
                                         report=F))
    karma<-summary(karma1)
    uwarma1<-(UWARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,rho=quant,
                             report=F))
    uwarma<-summary(uwarma1)
    barmax<-summary(BARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,
                                xreg = X,
                                report=F))
    karmax1<-suppressWarnings(KARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,
                          xreg = X,rho=quant,
                          control = list(method="Nelder-Mead",stopcr=1e-2),
                          report=F))
    karmax<-summary(karmax1)
    uwarmax1<-suppressWarnings(UWARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,
                                  rho=quant, xreg = X,
                                  report=F))
    uwarmax<-summary(uwarmax1)
    if(karma1$convergence==1 || is.nan(karma$aic)==1) karma$aic=0
    if(karmax1$convergence==1 || is.nan(karmax$aic)==1) karmax$aic=0
    if(uwarmax1$convergence==1 || is.nan(uwarmax$aic)==1) uwarmax$aic=0
    #   print(c(karma1$convergence,karma$aic))
    order[cont,]<-c(i,j,barma$aic,karma$aic,uwarma$aic,
                    barmax$aic,karmax$aic,uwarmax$aic)
    cont<-cont+1
  }
}
order<-order[-1,]

# ubxiix_best <- best.ubxii(hum_train, 
#                          sf = c(start = c(2020,356*24), frequency = 24*366),
#                          pmax = 3, qmax = 3,h=0,
#                          nbest = 1,X=X,X_hat = 0)
# 
# ubxii_best <- best.ubxii(hum_train, 
#                         sf = c(start = c(2020,356*24), frequency = 24*366),
#                         pmax = 3, qmax = 3,
#                         nbest = 1)

# selecting the order of each class
orbarma<-order[which(order[,3]==min(order[,3])),c(1:2)]
orkarma<-order[which(order[,4]==min(order[,4])),c(1:2)]
oruwarma<-order[which(order[,5]==min(order[,5])),c(1:2)]
orbxiiarmax<-c(3,0)
orbarmax<-order[which(order[,6]==min(order[,6])),c(1:2)]
orkarmax<-order[which(order[,7]==min(order[,7])),c(1:2)]
oruwarmax<-order[which(order[,8]==min(order[,8])),c(1:2)]
orbxiiarma<-c(1,2)

names_rows<-c("BARMAX","KARMAX","UWARMAX",
              "UBXII-ARMAX","ARIMAX",
              "BARMA","KARMA","UWARMA",
              "UBXII-ARMAX","ARIMA")

barmax<-BARFIMA.fit(hum_train,p=orbarmax[1],d=F,q=orbarmax[2],
                    xreg=X,info=T,report=F)
karmax<-KARFIMA.fit(hum_train,p=orkarmax[1],d=F,q=orkarmax[2],rho=quant,
                    control = list(method="Nelder-Mead",stopcr=1e-2),
                    xreg=X,info=T,report=F)
uwarmax<-UWARFIMA.fit(hum_train,p=oruwarmax[1],d=F,q=oruwarmax[2],rho=quant,
                      xreg=X,info=T,report=F)
ubxiiarmax<-ubxiiarma.fit(ts(hum_train),ar=1:3,ma=NA,
                         X=X,X_hat = 0,h=0)
barma<-BARFIMA.fit(hum_train,p=orbarma[1],d=F,q=orbarma[2],
                   info=T,report=F)
karma<-KARFIMA.fit(hum_train,p=orkarma[1],d=F,q=orkarma[2],rho=quant,
                   control = list(method="Nelder-Mead",stopcr=1e-2),
                   info=T,report=F)
uwarma<-UWARFIMA.fit(hum_train,p=oruwarma[1],d=F,q=oruwarma[2],rho=quant,
                     info=T,report=F)
ubxiiarma<-ubxiiarma.fit(ts(hum_train),ar=1,ma=1:2)

barmax_coeff<-t(summary(barmax)$coefficients[,c(1,4)])
karmax_coeff<-t(summary(karmax)$coefficients[,c(1,4)])
uwarmax_coeff<-t(summary(uwarmax)$coefficients[,c(1,4)])
ubxiiarmax_coeff<-t(ubxiiarmax$model[,c(1,4)])
arimax_coeff<-t(lmtest::coeftest(a02)[,c(1,4)])
barma_coeff<-t(summary(barma)$coefficients[,c(1,4)])
karma_coeff<-t(summary(karma)$coefficients[,c(1,4)])
uwarma_coeff<-t(summary(uwarma)$coefficients[,c(1,4)])
arima_coeff<-t(lmtest::coeftest(a01)[,c(1,4)])
ubxiiarma_coeff<-t(ubxiiarma$model[,c(1,4)])


xtable::xtable(round(barmax_coeff,4))
xtable::xtable(round(karmax_coeff,4))
xtable::xtable(round(uwarmax_coeff,4))
xtable::xtable(round(ubxiiarmax_coeff,4))
xtable::xtable(round(arimax_coeff[,c(5,1:4)],4))
xtable::xtable(round(barma_coeff,4))
xtable::xtable(round(karma_coeff,4))
xtable::xtable(uwarma_coeff, digits = 4)
xtable::xtable(round(ubxiiarma_coeff,4))
xtable::xtable(round(arima_coeff,4))

# Box.test BARMA
Box.test(barmax$residuals, lag=20, fitdf = sum(orbarmax))
# Box.test KARMA
Box.test(karmax$residuals, lag=20, fitdf = sum(orkarmax))
# Box.test BARMA
Box.test(uwarmax$residuals, lag=20, fitdf = sum(oruwarmax))
# Box.test BARMA
Box.test(ubxiiarmax$residuals, lag=20, fitdf = sum(oruwarmax))
# Box.test ARIMAX
Box.test(a02$residuals, lag=20, fitdf = 4)
# Box.test BARMA
Box.test(barma$residuals, lag=20, fitdf = sum(orbarma))
# Box.test KARMA
Box.test(karma$residuals, lag=20, fitdf = sum(orkarma))
# Box.test BARMA
Box.test(uwarma$residuals, lag=20, fitdf = sum(oruwarma))
# Box.test UBXII-ARMA
Box.test(ubxiiarma$residuals, lag=20, fitdf = 3)
# Box.test ARIMA
Box.test(a01$residuals, lag=20, fitdf = 4)

# w1=4.5
# h11=4
# setEPS()
# postscript("acf_beta.eps",width = w1, height = h11,family = "Times")
# par(mar = c(4.1, 4.1, 1, 1))
# acf(barmax$residuals, main ="")
# dev.off()
# 
# postscript("pacf_beta.eps",width = w1, height = h11,family = "Times")
# par(mar = c(4.1, 4.1, 1, 1))
# pacf(barmax$residuals, main ="")
# dev.off()
# 
# postscript("acf_KW.eps",width = w1, height = h11,family = "Times")
# par(mar = c(4.1, 4.1, 1, 1))
# acf(karmax$residuals, main ="")
# dev.off()
# 
# postscript("pacf_KW.eps",width = w1, height = h11,family = "Times")
# par(mar = c(4.1, 4.1, 1, 1))
# pacf(karmax$residuals, main ="")
# dev.off()
# 
# 
# postscript("acf_UW.eps",width = w1, height = h11,family = "Times")
# par(mar = c(4.1, 4.1, 1, 1))
# acf(uwarmax$residuals, main ="")
# dev.off()
# 
# postscript("pacf_UW.eps",width = w1, height = h11,family = "Times")
# par(mar = c(4.1, 4.1, 1, 1))
# pacf(uwarmax$residuals, main ="")
# dev.off()
# 
# postscript("acf_ubxii.eps",width = w1, height = h11,family = "Times")
# par(mar = c(4.1, 4.1, 1, 1))
# acf(ubxiiarmax$residuals, main ="")
# dev.off()
# 
# postscript("pacf_ubxii.eps",width = w1, height = h11,family = "Times")
# par(mar = c(4.1, 4.1, 1, 1))
# pacf(ubxiiarmax$residuals, main ="")
# dev.off()
# 
# postscript("acf_arima.eps",width = w1, height = h11,family = "Times")
# par(mar = c(4.1, 4.1, 1, 1))
# acf(a02$residuals, main ="")
# dev.off()
# 
# postscript("pacf_arima.eps",width = w1, height = h11,family = "Times")
# par(mar = c(4.1, 4.1, 1, 1))
# pacf(a02$residuals, main ="")
# dev.off()



barmax_out<-BARFIMA.extract(yt=hum, xreg = X0,  
                            coefs = list(alpha = barmax$coefficients[1], 
                                         beta = barmax$coefficients[2:(nX+1)],
                                         phi= barmax$coefficients[(nX+2):(orbarmax[1]+nX+1)], 
                                         theta = if(orbarmax[2]==0) {NULL} else{
                                           barmax$coefficients[(orbarmax[1]+(nX+2)):(orbarmax[1]+nX+1+orbarmax[2])]},
                                         nu = barmax$coefficients[(orbarmax[1]+(nX+2)+orbarmax[2])])
)
karmax_out<-KARFIMA.extract(yt=hum,xreg = X0,rho=quant,  
                            coefs = list(alpha = karmax$coefficients[1], 
                                         beta = karmax$coefficients[2:(nX+1)],
                                         phi= karmax$coefficients[(nX+2):(orkarmax[1]+nX+1)], 
                                         theta = if(orkarmax[2]==0) {NULL} else{
                                           karmax$coefficients[(orkarmax[1]+nX+2):(orkarmax[1]+nX+1+orkarmax[2])]},
                                         nu = karmax$coefficients[(orkarmax[1]+nX+2+orkarmax[2])])
)
uwarmax_out<-UWARFIMA.extract(yt=hum,xreg = X0,rho=quant,  
                              coefs = list(alpha = uwarmax$coefficients[1],
                                           beta = uwarmax$coefficients[2:(nX+1)],
                                           phi= uwarmax$coefficients[(nX+2):(oruwarmax[1]+nX+1)],
                                           theta = if(oruwarmax[2]==0) {NULL} else{
                                             uwarmax$coefficients[(oruwarmax[1]+nX+2):(oruwarmax[1]+nX+1+oruwarmax[2])]},
                                           nu = uwarmax$coefficients[(oruwarmax[1]+nX+2+oruwarmax[2])])
)
ubxiiarmax_out<-UWARFIMA.extract(yt=hum,xreg = X0,rho=quant,  
                              coefs = list(alpha = ubxiiarmax$alpha,
                                           beta = ubxiiarmax$beta,
                                           phi= ubxiiarmax$phi,
                                           theta = if(orbxiiarmax[2]==0) {NULL} else{
                                             ubxiiarmax$theta},
                                           nu = ubxiiarmax$c_par)
)
barma_out<-BARFIMA.extract(yt=hum,
                           coefs = list(alpha = barma$coefficients[1], 
                                        phi= barma$coefficients[2:(orbarma[1]+1)], 
                                        theta = if(orbarma[2]==0) {NULL} else{
                                          barma$coefficients[(orbarma[1]+2):(orbarma[1]+1+orbarma[2])]},
                                        nu = barma$coefficients[(orbarma[1]+2+orbarma[2])])
)

karma_out<-KARFIMA.extract(yt=hum,
                           coefs = list(alpha = karma$coefficients[1], 
                                        phi= karma$coefficients[2:(orkarma[1]+1)],
                                        theta = if(orkarma[2]==0) {NULL} else{
                                          karma$coefficients[(orkarma[1]+2):(orkarma[1]+1+orkarma[2])]},
                                        nu = karma$coefficients[(orkarma[1]+2+orkarma[2])])
)
uwarma_out<-UWARFIMA.extract(yt=hum,
                             coefs = list(alpha = uwarma$coefficients[1], 
                                          phi= uwarma$coefficients[2:(oruwarma[1]+1)],
                                          theta = if(oruwarma[2]==0) {NULL} else{
                                            uwarma$coefficients[(oruwarma[1]+2):(oruwarma[1]+1+oruwarma[2])]},
                                          nu = uwarma$coefficients[(oruwarma[1]+2+oruwarma[2])])
)

ubxiiarma_out<-UWARFIMA.extract(yt=hum,
                             coefs = list(alpha = ubxiiarma$alpha, 
                                          phi= ubxiiarma$phi,
                                          theta = ubxiiarma$theta,
                                          nu = ubxiiarma$c_par)
)

results_outsample<-rbind(
  forecast::accuracy(barmax_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(karmax_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(uwarmax_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(ubxiiarmax_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(new2$fitted, hum_test),
  forecast::accuracy(barma_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(karma_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(uwarma_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(ubxiiarma_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(new1$fitted, hum_test)
)[,c(3,2,5)]

row.names(results_outsample)<-
  names_rows

# xtable::xtable((results_outsample),digits=4)

# round(results_insample[,],5)
round(results_outsample[,],5)

df_percent_diff<-sweep(sweep(results_outsample, 2, results_outsample[3, ], FUN = "-"),  
                       2, results_outsample[3, ], FUN = "/")*100

df<-data.frame(values=as.vector(df_percent_diff),
               model=rep(names_rows,3),
               measure=rep(c("MAE","MAPE","RMSE"),10)
)
# w1=4.5
# h11=4
# setEPS()
# postscript("comparision1.eps",width = w1, height = h11,family = "Times")
ggplot(df,
       aes(y=values,x=measure,
           fill=factor(model))
)+
  geom_bar(stat="identity",position = "dodge", width = .9) +
  labs(fill="",y="Percentage differences",x="") +
  scale_fill_manual(values=gray.colors(9))+
  ylim(-1,480)+
  geom_text(aes(x=measure,y=values,
                label=round(values,2)),
            fontface="bold",
            hjust= 0,
            vjust = ifelse(df$values >= 0, -0.3, -.9),  # Ajusta a posição vertical dos números para ficarem acima da barra
            angle=55,
            position = position_dodge(width = .9),
            color="black",
            size=2)+
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold",size=8),
        plot.title = element_text(face="bold",size=8),
        legend.text = element_text(face="bold",size=8),
        legend.key.size = unit(0.3, "cm"),  # Ajusta o tamanho dos quadrados na legenda
        legend.spacing.x = unit(0.3, 'cm'),  # Espaçamento entre os itens da legenda
        legend.margin = margin(t = -13, unit = "pt"),
        axis.title.y = element_text(face="bold", color="black",
                                    size=8),
        axis.title.x = element_text(face="bold", color="black",
                                    size=8),
        axis.text.x = element_text(face="bold", color="black",
                                   size=8),
        axis.text.y = element_text(face="bold", color="black",
                                   size=8),
        panel.background = element_rect(fill = "white", colour = "white"))

# dev.off()

# hum1=xts::xts(hum_test, order.by=data$timestamp[(n+1):(dim(data)[1])])
# hum2<-xts::xts(new2$fitted, order.by=data$timestamp[(n+1):(dim(data)[1])])#
# hum3<-xts::xts(uwarmax_out$mut[(n+1):(dim(data)[1])], order.by=data$timestamp[(n+1):(dim(data)[1])])
# w1=4.5
# h11=2.5
# setEPS()
# postscript("comparision2.eps",width = w1, height = h11,family = "Times")
# plot(hum1,
#      main="", yaxis.right=FALSE, grid.col = "white",
#      format.labels="%b-%Y", main.timespan = FALSE,
#      lwd=.5,
#      ylim=c(min(hum1,hum2,hum3),
#             max(hum1,hum2,hum3)+.1)
# )
# lines(hum2,lty=2,lwd=1,col=2)
# lines(hum3,lwd=.5,col=4)
# xts::addLegend("topleft",
#                legend.names =c("Original data","ARIMAX","UWARMAX"),
#                col=c(1,2,4), cex=.8,lty=c(1,2,1),
#                lwd=c(.5,1,.5),
#                ncol=3
# )
# dev.off()
