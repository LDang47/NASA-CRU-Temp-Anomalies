rm(list = ls())
x <- c('dplyr','lubridate','plotly','forecast','tseries','MLmetrics') 
lapply(x,require,character.only=TRUE)

#Monthly Temperature - Global (NH+SH)/2
#Time series are presented as temperature anomalies (deg C) relative to 1961-1990.


df <- read.table(url("https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.6.0.0.monthly_ns_avg.txt"), 
                           header=FALSE,fill = T,row.names = NULL )
names(df) <- c("Date","value",
                         "LB_95_CI_bias_uncertainty","UB_95_CI_bias_uncertainty",
                         "LB_95_CI_measurement_sampling_uncertainty","UB_95_CI_measurement_sampling_uncertainty",
                         "LB_95_CI_coverage_uncertainty","UB_95_CI_coverage_uncertainty",
                         "LB_95_CI_measurement_sampling_bias_uncertainty","UB_95_CI_measurement_sampling_bias_uncertainty",
                         "LB_95_CI_measurement_sampling_bias_coverage_uncertainty","UB_95_CI_measurement_sampling_bias_coverage_uncertainty")



# Dropping not required columns
drop_col <- c("LB_95_CI_bias_uncertainty","UB_95_CI_bias_uncertainty",
              "LB_95_CI_measurement_sampling_uncertainty","UB_95_CI_measurement_sampling_uncertainty",
              "LB_95_CI_coverage_uncertainty","UB_95_CI_coverage_uncertainty",
              "LB_95_CI_measurement_sampling_bias_uncertainty","UB_95_CI_measurement_sampling_bias_uncertainty",
              "LB_95_CI_measurement_sampling_bias_coverage_uncertainty","UB_95_CI_measurement_sampling_bias_coverage_uncertainty")

df <- df %>% select(-c(drop_col))

df

df_monthly <- ts(df$value,start=1850, frequency=12)
plot(df_monthly, xlab="Year-Month", ylab="Temperature Celsius")

CRU_dts <- decompose(df_monthly)
plot(CRU_dts)



### Check for model accuracy

CRU_msts <- msts (df$value, seasonal.periods=c(12,120)) # #Yearly and Decade
CRU_msts_120_240 <- msts (df$value, seasonal.periods=c(12,120,240)) # # Yearly Decade and 2*Decade
CRU_msts_240 <- msts (df$value, seasonal.periods=c(12,240)) # # Yearly and 2*Decade
class(df$value)

accuracy_df <- data.frame(model=c(""),runid=c(""),metric_mape=c(as.double("")),metric_mae=c(as.double("")))
# THE LOOP IS FOR ONLY MS_TS OBJECT MODELS
for (i in 1:10)
{ 
  nTest <- 12*i  
  nTrain <- length(CRU_msts)- nTest -1
  train <- window(CRU_msts, start=1, end=1+(nTrain)/(120))
  test <- window(CRU_msts, start=1+(nTrain+1)/(120), end=1+(nTrain+12)/(120))
  
  cat("PRE RUN COMMENT----------------------------------
      Data Partition",i,"
      Training Set includes",nTrain," time periods. Observations 1 to", nTrain, "
      Test Set includes 12 time periods. Observations", nTrain+1, "to", nTrain+12,"
      ")
  # Doing Fourier ARIMA
  print("Running Fourier ARIMA")
  arima_fourier <- auto.arima(train, seasonal=FALSE,
                    xreg=fourier(train, K=c(5,5)))
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("ARIMA_FOURIER_12_120"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(forecast(arima_fourier,xreg=fourier(test, K=c(5,5)),h=length(test))[['mean']],test)
                                              )),
                                            metric_mae=c(paste(
                                              MAE(forecast(arima_fourier,xreg=fourier(test, K=c(5,5)),h=length(test))[['mean']],test)
                                            ))))
  print("END Fourier ARIMA")
  # Doing STL + ETS
  print("Running STL ETS")
  stl_forecast <- stlf(train,h=length(test))[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("STL_ETS_12_120"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(stl_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(stl_forecast,test)
                                            ))
                                            ))

  # Doing STL ARIMA
  print("Running STL ARIMA")
  stl_forecast <- stlf(train,h=length(test),method='arima')[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("STL_ARIMA_12_120"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(stl_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(stl_forecast,test)
                                            ))
                                            ))
  # TBATS
  print("Running TBATS")
  tbats_forecast <- forecast(tbats(train), h=length(test))[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("TBATS_12_120"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(tbats_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(tbats_forecast,test)
                                            ))
                                            ))
  

  }

accuracy_df <- accuracy_df[-1,]
accuracy_df$metric_mape <- as.double(accuracy_df$metric_mape)
accuracy_df$metric_mae <- as.double(accuracy_df$metric_mae)
accuracy_mean <- accuracy_df %>% dplyr::group_by(model) %>% summarise(mape_mean=mean(metric_mape),mae_mean=mean(metric_mae))



for (i in 1:20)
{ 
  print(Sys.time())
  nTest <- 12*i  
  nTrain <- length(nasa_msts_120_240)- nTest -1
  train <- window(nasa_msts_120_240, start=1, end=1+(nTrain)/(240))
  test <- window(nasa_msts_120_240, start=1+(nTrain+1)/(240), end=1+(nTrain+12)/(240))
  
  
  cat("PRE RUN COMMENT----------------------------------
      Data Partition",i,"
      Training Set includes",nTrain," time periods. Observations 1 to", nTrain, "
      Test Set includes 12 time periods. Observations", nTrain+1, "to", nTrain+12,"
      ")
  
  # Doing Fourier ARIMA
  print("Running Fourier ARIMA")
  arima_fourier <- auto.arima(train, seasonal=FALSE,
                              xreg=fourier(train, K=c(5,5,5)))
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("ARIMA_FOURIER_12_120_240"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(forecast(arima_fourier,xreg=fourier(test, K=c(5,5,5)),h=length(test))[['mean']],test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(forecast(arima_fourier,xreg=fourier(test, K=c(5,5,5)),h=length(test))[['mean']],test)
                                            ))
                                            ))
  print("END Fourier ARIMA")
  # Doing STL + ETS
  print("Running STL ETS")
  stl_forecast <- stlf(train,h=length(test))[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("STL_ETS_12_120_240"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(stl_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(stl_forecast,test)
                                            ))
                                            ))
  
  # Doing STL ARIMA
  print("Running STL ARIMA")
  stl_forecast <- stlf(train,h=length(test),method='arima')[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("STL_ARIMA_12_120_240"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(stl_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(stl_forecast,test)
                                            ))
                                            ))
  # TBATS
  print("Running TBATS")
  tbats_forecast <- forecast(tbats(train), h=length(test))[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("TBATS_12_120_240"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(tbats_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(tbats_forecast,test)
                                            ))
                                            ))
  
}


accuracy_df$metric_mape <- as.double(accuracy_df$metric_mape)
accuracy_df$metric_mae <- as.double(accuracy_df$metric_mae)
accuracy_mean <- accuracy_df %>% dplyr::group_by(model) %>% summarise(mape_mean=mean(metric_mape),mae_mean=mean(metric_mae))


## For series 12, 240
for (i in 1:20)
{ 
  print(Sys.time())
  nTest <- 12*i  
  nTrain <- length(nasa_msts_240)- nTest -1
  train <- window(nasa_msts_240, start=1, end=1+(nTrain)/(240))
  test <- window(nasa_msts_240, start=1+(nTrain+1)/(240), end=1+(nTrain+12)/(240))
  
  
  cat("PRE RUN COMMENT----------------------------------
      Data Partition",i,"
      Training Set includes",nTrain," time periods. Observations 1 to", nTrain, "
      Test Set includes 12 time periods. Observations", nTrain+1, "to", nTrain+12,"
      ")
  
  # Doing Fourier ARIMA
  print("Running Fourier ARIMA")
  arima_fourier <- auto.arima(train, seasonal=FALSE,
                              xreg=fourier(train, K=c(5,5)))
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("ARIMA_FOURIER_12_240"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(forecast(arima_fourier,xreg=fourier(test, K=c(5,5)),h=length(test))[['mean']],test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(forecast(arima_fourier,xreg=fourier(test, K=c(5,5)),h=length(test))[['mean']],test)
                                            ))
                                            ))
  print("END Fourier ARIMA")
  # Doing STL + ETS
  print("Running STL ETS")
  stl_forecast <- stlf(train,h=length(test))[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("STL_ETS_12_240"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(stl_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(stl_forecast,test)
                                            ))
                                            ))
  
  # Doing STL ARIMA
  print("Running STL ARIMA")
  stl_forecast <- stlf(train,h=length(test),method='arima')[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("STL_ARIMA_12_240"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(stl_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(stl_forecast,test)
                                            ))
                                            ))
  
  # TBATS
  print("Running TBATS")
  tbats_forecast <- forecast(tbats(train), h=length(test))[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("TBATS_12_240"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(tbats_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(tbats_forecast,test)
                                            ))
                                            ))
}


accuracy_df$metric_mape <- as.double(accuracy_df$metric_mape)
accuracy_df$metric_mae <- as.double(accuracy_df$metric_mae)
accuracy_mean <- accuracy_df %>% dplyr::group_by(model) %>% summarise(mape_mean=mean(metric_mape),mae_mean=mean(metric_mae))


# TS of 12 seasonality
for (i in 1:10)
{ 
  print(Sys.time())
  train <- ts(df$value,start = c(1880,1),end=c((2019-i),12),frequency = 12)
  test <- ts(df$value[length(train)+1:12],start = c((2019-i+1),1),frequency = 12)
  
  cat("PRE RUN COMMENT----------------------------------
      Data Partition",i,"
      Training Set includes",length(train)," time periods. Observations 1 to", length(train), "
      Test Set includes 12 time periods. Observations", length(train)+1, "to", length(train)+12,"
      ")
  # Doing Fourier ARIMA
  print("Running Fourier ARIMA")
  arima_fourier <- auto.arima(train, seasonal=FALSE,
                              xreg=fourier(train, K=c(5)))
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("ARIMA_FOURIER_12"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(forecast(arima_fourier,xreg=fourier(test, K=c(5)),h=length(test))[['mean']],test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(forecast(arima_fourier,xreg=fourier(test, K=c(5)),h=length(test))[['mean']],test)
                                            ))
                                            ))
  print("END Fourier ARIMA")
  
  # Doing Plain Vanilla ARIMA
  print("Running Plain Vanilla ARIMA")
  arima_d_1 <- auto.arima(train, D=1)
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("ARIMA_1_D_12"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(forecast(arima_d_1,h=length(test))[['mean']],test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(forecast(arima_d_1,h=length(test))[['mean']],test)
                                            ))
  ))
  print("END Plain Vanilla ARIMA")
  # Doing STL + ETS
  print("Running STL ETS")
  stl_forecast <- stlf(train,h=length(test))[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("STL_ETS_12"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(stl_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(stl_forecast,test)
                                            ))
                                            ))
  
  # Doing STL ARIMA
  print("Running STL ARIMA")
  stl_forecast <- stlf(train,h=length(test),method='arima')[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("STL_ARIMA_12"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(stl_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(stl_forecast,test)
                                            ))
                                            ))
  # TBATS
  print("Running TBATS")
  tbats_forecast <- forecast(tbats(train), h=length(test))[['mean']]
  accuracy_df<-rbind(accuracy_df,data.frame(model=c("TBATS_12"),
                                            runid=c(paste(i)),
                                            metric_mape=c(paste(
                                              MAPE(tbats_forecast,test)
                                            )),
                                            metric_mae=c(paste(
                                              MAE(tbats_forecast,test)
                                            ))
                                            ))
  
  
}


accuracy_df$metric_mape <- as.double(accuracy_df$metric_mape)
accuracy_df$metric_mae <- as.double(accuracy_df$metric_mae)
accuracy_mean <- accuracy_df %>% dplyr::group_by(model) %>% summarise(mape_mean=mean(metric_mape),mae_mean=mean(metric_mae))


###################################################################
#######ANSWERING QUESTIONS#########################################


# Create multiple seasonalities: one for 12 months and 10 years and another for 12 months and 20 years

######Analysis with Anomalies


df$value
df_monthly <- ts(df$value,start=1850, frequency=12)
plot(df_monthly, xlab="Year-Month", ylab="Anomalies",main="CRU Monthly Data")


CRU_msts <- msts (df$value, seasonal.periods=c(12,120)) # #Yearly and Decade
CRU_msts_120_240 <- msts (df$value, seasonal.periods=c(12,120,240)) # # Yearly Decade and 2*Decade
CRU_msts_240 <- msts (df$value, seasonal.periods=c(12,240)) # # Yearly and 2*Decade
class(df$value)


#Plot for the seasonalities

tbats1 <- tbats(CRU_msts)
plot(tbats1) #plot decomposition

plot(tbats(CRU_msts))

plot(mstl(CRU_msts))
plot(mstl(CRU_msts_120_240))
plot(mstl(CRU_msts_240))

plot(tbats(CRU_msts))
plot(tbats(CRU_msts_120_240))
plot(tbats(CRU_msts_240))


stlf(CRU_msts)


###Testing

df

kpss.test(df_monthly)
mean(df_monthly)
sd(df_monthly)
var(df_monthly)
Box.test(df_monthly,type="Ljung")


#Arima Fourier 12_120 model

arima_fourier <- auto.arima(CRU_msts, seasonal=FALSE,
                            xreg=fourier(CRU_msts, K=c(5,5)))

checkresiduals(arima_fourier)

#forecast the model for the next 972 months

arimaf_model<-forecast(arima_fourier,xreg=fourier(CRU_msts,K=c(5,5),h=972))
plot(arimaf_model)

#view the model

View(data.frame(forecast(arima_fourier,xreg=fourier(CRU_msts,K=c(5,5),h=972))))

#Send the model to excel 

write.csv(arimaf_model, file = "arima_prediction.csv") # export the selected model's predictions into a CSV file


#######ANALYSIS WITH TEMPERATURES#####################

df$valueT<-df$value+14

df_monthlyT <- ts(df$valueT,start=1850, frequency=12)
plot(df_monthlyT, xlab="Year-Month", ylab="Temperatures in Celsius")

View(df_monthlyT)

CRU_msts <- msts (df$valueT, seasonal.periods=c(12,120)) # #Yearly and Decade
CRU_msts_120_240 <- msts (df$valueT, seasonal.periods=c(12,120,240)) # # Yearly Decade and 2*Decade
CRU_msts_240 <- msts (df$valueT, seasonal.periods=c(12,240)) # # Yearly and 2*Decade
class(df$value)

#Plot for the seasonalities

plot(mstl(CRU_msts))
plot(mstl(CRU_msts_120_240))
plot(mstl(CRU_msts_240))
stlf(CRU_msts)

#Arima Fourier 12_120 model

arima_fourier <- auto.arima(CRU_msts, seasonal=FALSE,
                            xreg=fourier(CRU_msts, K=c(5,5)))

checkresiduals(arima_fourier)

#forecast the model for the next 972 months

arimaf_model<-forecast(arima_fourier,xreg=fourier(CRU_msts,K=c(5,5),h=972),level=c(80,90,95))
plot(arimaf_model)

#view the model

View(data.frame(forecast(arima_fourier,xreg=fourier(CRU_msts,K=c(5,5),h=972))))

#Send the model to excel 

write.csv(arimaf_model, file = "arima_prediction.csv") # export the selected model's predictions into a CSV file




##Bringing the NASA dataframe

dfnasa <- read.csv("https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv",skip = 1)
head(dfnasa)


# Dropping not required columns
drop_col <- c("J.D","D.N","DJF","MAM","JJA","SON")

dfnasa <- dfnasa %>% select(-c(drop_col))

dfnasa <- dfnasa %>% mutate_if(is.factor,as.character)
df_monthly_N <- reshape2::melt(dfnasa,id.var=c('Year'))              

df_monthly_N


## Dropping the value from 2020 April as it is populated as *****
df_monthly_N$Date <- mdy(paste(df_monthly_N$variable,'-01-',df_monthly_N$Year,sep=""))
df_monthly_N

df_monthly_N <- df_monthly_N %>% dplyr::filter(Date<'2020-01-01') %>% arrange(Date)


df_monthly_N$value <- as.numeric(df_monthly_N$value)
head(df_monthly_N)
plot(df_monthly_N$Date,df_monthly_N$value,type = "o")

# Lets test model
df_monthly_N$value<-df_monthly_N$value+14

nasa_ts <- ts(df_monthly_N$value,start = c(1880,1),frequency = 12)
plot(nasa_ts)
lines(CRU, col = "red")


#Taking the same intervals
df_monthly_N <- df_monthly_N %>% dplyr::filter(Date<'2020-01-01') %>% arrange(Date)

CRU <- ts(df$valueT,start=1850, frequency=12)

ts.plot(nasa_ts, CRU, gpars = list(col = c("blue", "red")),xlab="Date", ylab="Temperature in Celsius", lty=c(1:3))



mean(CRU)
mean(nasa_ts)


