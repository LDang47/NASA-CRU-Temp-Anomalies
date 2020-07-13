rm(list = ls())
x <- c('dplyr','lubridate','plotly','forecast','tseries','MLmetrics') 
lapply(x,require,character.only=TRUE)

df <- read.csv("https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv",skip = 1)
# Dropping not required columns
drop_col <- c("J.D","D.N","DJF","MAM","JJA","SON")

df <- df %>% select(-c(drop_col))

df <- df %>% mutate_if(is.factor,as.character)
df_monthly <- reshape2::melt(df,id.var=c('Year'))              

## Dropping the value from 2020 April as it is populated as *****
df_monthly$Date <- mdy(paste(df_monthly$variable,'-01-',df_monthly$Year,sep=""))


df_monthly <- df_monthly %>% dplyr::filter(Date<'2020-01-01') %>% arrange(Date)

df_monthly$value <- as.numeric(df_monthly$value)
plot(df_monthly$Date,df_monthly$value,type = "o")

# Reference Links:
# https://data.giss.nasa.gov/gistemp/faq/#q101
# Download the data ffrom the below link  ; ; Search for 'Global-mean monthly, seasonal, and annual means'
# https://data.giss.nasa.gov/gistemp/


nasa_ts <- ts(df_monthly$value,start = c(1880,1),frequency = 12)
autoplot(nasa_ts)
nasa_dts <- decompose(nasa_ts)
plot(nasa_dts)

dygraph(nasa_ts,main="NASA Global Anomalies Temp Time Series")





plot(stl(nasa_ts,s.window=12))
  
monthplot(nasa_ts)
forecast::seasonplot(nasa_ts)
plot(cycle(nasa_ts))  

plot(decompose(nasa_ts,type = c("multiplicative")))  
plot(stl(nasa_ts,s.window = "periodic"),main="NASA Climate")

nasa_ar_model <- auto.arima(nasa_ts,seasonal = TRUE)
nasa_ar_model <- auto.arima(nasa_msts,D=1)
plot(forecast(nasa_ar_model,960))  
ggtsdisplay(nasa_ts)
checkresidual(nasa_ts)


nasa_ar_model_2 <- auto.arima(nasa_ts,D=2)

plot(forecast(nasa_ar_model,960))  
ggtsdisplay(nasa_ts)
checkresidual(nasa_ts)

nasa_nnet <- nnetar(nasa_ts)

nnetforecast <- forecast(nasa_nnet, h = 1000,
                         PI = T)

plot(nnetforecast)
nasa_msts <- msts (df_monthly$value, seasonal.periods=c(12,120))
nasa_msts_120_240 <- msts (df_monthly$value, seasonal.periods=c(12,120,240))
nasa_msts_240 <- msts (df_monthly$value, seasonal.periods=c(12,240))

nasa_tbats <- tbats(nasa_msts)
nasa_tbats_12_24 <- tbats(nasa_msts_120_240)
plot(nasa_tbats_12_24) #plot decomposition
plot(nasa_tbats) #plot decomposition
plot(tbats(nasa_msts_240)) #plot decomposition

nasa_tbats_pred <- forecast(nasa_tbats, h=960, level=c(0.8, 0.95)) #predict 2 weeks out
plot(nasa_tbats_pred, xlab="Time", ylab="Predicted Anaomolies")



##############
nasalm_msts <- tslm(nasa_msts ~ trend + season) # Build a linear model for trend and seasonality
summary(nasalm_msts)

residarima1 <- auto.arima(nasalm_msts$residuals) # Build ARIMA on it's residuals
residarima1
residualsArimaForecast <- forecast(residarima1, h=960) #forecast from ARIMA
residualsF <- as.numeric(residualsArimaForecast$mean)

regressionForecast <- forecast(nasalm_msts,h=960) #forecast from lm
regressionF <- as.numeric(regressionForecast$mean)

forecastR <- regressionF+residualsF # Total prediction
print(forecastR)
for (i in 1:960){points(i+1680,forecastR[i],col="red",pch=19, cex=0.5)}

#compare with TBATS
plot(forecast(arima_fourier,xreg=fourier(nasa_msts, K=c(5,5), h=960)), xlab="Time", ylab="Predicted anamolies")
for (i in 1:960){points((i+1680+120)/(120),forecastR[i],col="red",pch=19, cex=0.5)}





###########################
# This from the FPP BOOK
###########################

#STL with multiple seasonal periods 

plot(mstl(nasa_msts))
plot(mstl(nasa_msts_120_240))
stlf(nasa_msts)
plot(stlf(nasa_msts,h=972))
plot(stlf(nasa_msts,h=972,method="arima"))
plot(stlf(nasa_msts,h=960))
plot(stlf(nasa_msts_120_240,h=960))

##Dynamic harmonic regression with multiple seasonal periods
arima_fourier <- auto.arima(nasa_msts, seasonal=FALSE,
                  xreg=fourier(nasa_msts, K=c(5,5)))

plot(forecast(arima_fourier,xreg=fourier(nasa_msts, K=c(5,5), h=972)))

checkresiduals(fit)

fit_man <- Arima(nasa_msts, order=c(3,1,3),# lambda=0
                  xreg=fourier(nasa_msts, K=c(5,5)))

plot(forecast(fit_man,xreg=fourier(nasa_msts, K=c(5,5), h=960)))

checkresiduals(fit)


ks_arima_fourier_1 <- Arima(ks_msts, order=c(3,1,3), 
                            xreg=fourier(ks_msts, K=c(5,5)))



ks_arimafou_pred_1 <- forecast(ks_arima_fourier_1, level=c(0.8,0.9),xreg=fourier(ks_msts, K=c(5,5), h=972))

# TBATS models
nasa_ms_tbats <- tbats(nasa_msts)
plot(forecast(nasa_ms_tbats, h=972))

# Neuraal network with multiple seasonality
nasa_ms_nnet <- nnetar(nasa_msts)
autoplot(forecast(nasa_ms_nnet,PI=T,h=960))


simple_arima <- autoplot(forecast(nasa_ar_model,h=960))

##########################################
## FORECAST COMBINATION
##########################################
autoplot(nasa_msts) +
autolayer(forecast(nasa_ms_tbats,h=960,level = c(0)),series="MSTS TBATS") +
autolayer(forecast(fit,xreg=fourier(nasa_msts, K=c(5,5), h=960), level = c(0)),series = "FOURIER ARIMA") 
# autolayer(simple_arima[["mean"]], series="TS ARIMA")+
#   scale_x_date(date_breaks = "months" , date_labels = "%b-%y")
#   xlab("Year") + ylab("Temp")
#   ggtitle("TEMperature anamolies")


