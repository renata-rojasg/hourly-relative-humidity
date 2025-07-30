######################
## Data preparation ##
######################
library(tidyverse)
data <- readr::read_delim("https://raw.githubusercontent.com/emanueleg/lora-rssi/master/vineyard-2021_data/combined_hourly_data.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE) |> 
  mutate(timestamp=as.POSIXct(timestamp, tz="GMT",
                              origin="1970-01-01 00:00:00"),
         hum=hum/100)
#
n<-round(dim(data)[1]*.8)-149
data<-data[830:1870,]

#########################
## Descriptive analysis ##
#########################
length(table(strftime(as.POSIXct(data$timestamp, tz="GMT",
                                 origin="1970-01-01 00:00:00"),
                      format="%Y%m%d", tz="GMT"))) # number of days

summary1<-apply(data[,c(3,2)],2,summary)
descriptive<-(rbind(summary1[c(1,3,4,6),],
                    SD=apply(data[,c(3,2)],2,sd),
                    CV=apply(data[,c(3,2)],2,sd)/
                      abs(apply(data[,c(3,2)],2,mean))*100,
                    Skew=apply(data[,c(3,2)],2,
                               moments::skewness)
                    ))
# Table 2
print(t(descriptive),digits=2)

###########
## plots ##
###########
# RH plots
RH <- xts::xts(data$hum, order.by=data$timestamp)

#Fig 1a
RHplot<-{plot(RH, main="",
              yaxis.right=FALSE, grid.col = "white",
              format.labels="%b-%Y", main.timespan = FALSE,
              lwd=0.5,cex.lab=1.3,cex.axis=1.3,
              ylim=c(0.4,.9))
}

RHplot

#Fig 1b
hist(RH,main="",xlab = "",freq = T)

# Temp plots
Temp <- xts::xts(data$temp, order.by=data$timestamp)

Temp_plot<-{plot(Temp, main="",
              yaxis.right=FALSE, grid.col = "white",
              format.labels="%b-%Y", main.timespan = FALSE,
              lwd=0.5,cex.lab=1.3,cex.axis=1.3,
              ylim=c(-5,10))
}

#Fig 2a
Temp_plot

#Fig 2
hist(Temp,main="",xlab = "",freq = T)
