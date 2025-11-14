# rm(list = ls())
library(forecast)
library(dplyr)
library(BTSR)
# source("Section3_2_test05_complete.R")
# source("Section3_2_test05_dummy.R")
# source("Section3_2_test05_autumn.R")
source("Section_3_1 - completo.R")
#############################
## Missing data experiment ##
#############################
# fitting the algorithms
steps_values<-round(length(hum_test)*c(.01,.05,.1))
steps<-length(steps_values)
meanRMSE<-meanMAE<-meanMAPE<-
  rRMSE<-rMAE<-rMAPE<-matrix(NA,steps,3)

nx=c(2)

for(j in 1:steps){
  n_ahead<-steps_values[j]
  hat1<-hat2<-hat3<-matrix(NA,length(hum_test)-n_ahead+1,3)
  colnames(hat1)<-colnames(hat2)<-colnames(hat3)<-c("RMSE","MAE","MAPE")
  for(k in 0:(length(hum_test)-n_ahead)){
    # updating the data
    hum_test1<-hum[1:(n+k)]
    n_new<-length(hum_test1)
    X1<-as.matrix(data[1:(n+k),nx])
    Xtest1<-as.matrix(data[(n+k+1):(n+k+n_ahead),nx])
    # last n_ahead observations for the MA computation
    hum_prev<-hum_test1[(n_new-n_ahead+1):(n_new)] 
    # original values 
    hum_true<-hum[(n_new+1):(n_new+n_ahead)]
    # moving average
    hum_hat1<-rep(mean(hum_prev),n_ahead)
    # ARIMAX
    new02<-Arima(hum_test1,model=a02,xreg = X1)
    hum_hat2<-predict(new02,
                      n.ahead = n_ahead,
                      newxreg=Xtest1)$pred
    # UWARIMAX
    hum_hat3<- UWARFIMA.extract(yt=hum_test1,
                                xreg = X1,rho=quant,
                                xnew = Xtest1,
                                nnew = n_ahead,
                                coefs = list(alpha = uwarmax$coefficients[1],
                                             beta = uwarmax$coefficients[2:(nX+1)],
                                             phi= uwarmax$coefficients[(nX+2):(oruwarmax[1]+nX+1)],
                                             theta = if(oruwarmax[2]==0) {NULL} else{
                                               uwarmax$coefficients[(oruwarmax[1]+nX+2):(oruwarmax[1]+nX+1+oruwarmax[2])]},
                                             nu = uwarmax$coefficients[(oruwarmax[1]+nX+2+oruwarmax[2])]
                                )
    )$yt.new
    hat1[k+1,]<-forecast::accuracy(hum_hat1, hum_true)[c(2,3,5)]
    hat2[k+1,]<-forecast::accuracy(hum_hat2, hum_true)[c(2,3,5)]
    hat3[k+1,]<-forecast::accuracy(hum_hat3, hum_true)[c(2,3,5)]
  }
  # accuracy measures
  RMSE<-cbind(hat1[,1],hat2[,1],hat3[,1])
  MAE <-cbind(hat1[,2],hat2[,2],hat3[,2])
  MAPE<-cbind(hat1[,3],hat2[,3],hat3[,3])
  
  colnames(meanRMSE)<-colnames(meanMAE)<-
    colnames(meanMAPE)<-colnames(rRMSE)<-
    colnames(rMAE)<-colnames(rMAPE)<-
    c("MA-imputation", "ARIMAX imputation", "UWARIMAX imputation")
  rownames(meanRMSE)<-rownames(meanMAE)<-
    rownames(meanMAPE)<-rownames(rRMSE)<-
    rownames(rMAE)<-rownames(rMAPE)<-
    c("1%","5%","10%")
  meanRMSE[j,]<-apply(RMSE,2,mean)
  meanMAE[j,]<-apply(MAE,2,mean)
  meanMAPE[j,]<-apply(MAPE,2,mean)
  rRMSE[j,]<-apply(apply(RMSE, 1, rank)==1,1,sum)
  rMAE[j,]<-apply(apply(MAE, 1, rank)==1,1,sum)
  rMAPE[j,]<-apply(apply(MAPE, 1, rank)==1,1,sum)
}



RMSE_results<-data.frame(
  model=c(rep("MA",3),
          rep("ARIMAX",3),
          rep("UWARMAX",3)),
  gap=rep(c("1%","5%","10%"),3),
  value=as.vector(meanRMSE)
)

RMSE_results$model <- factor(RMSE_results$model, levels = c("MA", "ARIMAX", "UWARMAX"))


MAPE_results<-data.frame(
  model=c(rep("MA",3),
          rep("ARIMAX",3),
          rep("UWARMAX",3)),
  gap=rep(c("1%","5%","10%"),3),
  value=as.vector(meanMAPE)
)

MAPE_results$model <- factor(MAPE_results$model, levels = c("MA", "ARIMAX", "UWARMAX"))

MAE_results<-data.frame(
  model=c(rep("MA",3),
          rep("ARIMAX",3),
          rep("UWARMAX",3)),
  gap=rep(c("1%","5%","10%"),3),
  value=as.vector(meanMAE)
)

MAE_results$model <- factor(MAE_results$model, levels = c("MA", "ARIMAX", "UWARMAX"))

p1<-(ggplot(RMSE_results, aes(x=gap,y=value,fill=model))+
         geom_bar(stat="identity",position = "dodge")+
         labs(fill="",x="Gap sizes", y="RMSE")  +
         # ylim(0,4)+
         geom_text(aes(x=gap,y=value/2,
                       label=round(value,3)),
                   fontface="bold",
                   angle=90,
                   position = position_dodge(width = 1),
                   # color="black",
                   size=6)+
         scale_x_discrete(limits=c("1%","5%","10%"))+
         scale_fill_manual(values=c("grey90","grey65","grey45") )+
         theme(legend.position = "bottom",
               strip.text = element_text(face="bold",size=15),
               plot.title = element_text(face="bold",size=15),
               legend.text = element_text(face="bold",size=15),
               axis.title.y = element_text(face="bold", color="black",
                                           size=15),
               axis.title.x = element_text(face="bold", color="black",
                                           size=15),
               axis.text.x = element_text(face="bold", color="black",
                                          size=15),
               axis.text.y = element_text(face="bold", color="black",
                                          size=15),
               panel.background = element_rect(fill = "white", colour = "white")))

# p1

p2<-(ggplot(MAPE_results, aes(x=gap,y=value,fill=model))+
       geom_bar(stat="identity",position = "dodge")+
       labs(fill="",x="Gap sizes", y="MAPE")  +
       # ylim(0,4)+
       geom_text(aes(x=gap,y=value/2,
                     label=round(value,3)),
                 fontface="bold",
                 angle=90,
                 position = position_dodge(width = 1),
                 # color="black",
                 size=6)+
       scale_x_discrete(limits=c("1%","5%","10%"))+
       scale_fill_manual(values=c("grey90","grey65","grey45") )+
       theme(legend.position = "bottom",
             strip.text = element_text(face="bold",size=15),
             plot.title = element_text(face="bold",size=15),
             legend.text = element_text(face="bold",size=15),
             axis.title.y = element_text(face="bold", color="black",
                                         size=15),
             axis.title.x = element_text(face="bold", color="black",
                                         size=15),
             axis.text.x = element_text(face="bold", color="black",
                                        size=15),
             axis.text.y = element_text(face="bold", color="black",
                                        size=15),
             panel.background = element_rect(fill = "white", colour = "white")))

# p2

p3<-(ggplot(MAE_results, aes(x=gap,y=value,fill=model))+
       geom_bar(stat="identity",position = "dodge")+
       labs(fill="",x="Gap sizes", y="MAE")  +
       # ylim(0,4)+
       geom_text(aes(x=gap,y=value/2,
                     label=round(value,3)),
                 fontface="bold",
                 angle=90,
                 position = position_dodge(width = 1),
                 # color="black",
                 size=6)+
       scale_x_discrete(limits=c("1%","5%","10%"))+
       scale_fill_manual(values=c("grey90","grey65","grey45") )+
       theme(legend.position = "bottom",
             strip.text = element_text(face="bold",size=15),
             plot.title = element_text(face="bold",size=15),
             legend.text = element_text(face="bold",size=15),
             axis.title.y = element_text(face="bold", color="black",
                                         size=15),
             axis.title.x = element_text(face="bold", color="black",
                                         size=15),
             axis.text.x = element_text(face="bold", color="black",
                                        size=15),
             axis.text.y = element_text(face="bold", color="black",
                                        size=15),
             panel.background = element_rect(fill = "white", colour = "white")))

# p3

w1=4.5
h11=4
setEPS()
postscript("RMSE.eps",width = w1, height = h11,family = "Times")
p1
dev.off()
setEPS()
postscript("MAPE.eps",width = w1, height = h11,family = "Times")
p2
dev.off()

setEPS()
postscript("MAE.eps",width = w1, height = h11,family = "Times")
p3
dev.off()