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

# Lets test model

nasa_ts <- ts(df_monthly$value,start = c(1880,1),frequency = 12)
plot(nasa_ts)
nasa_dts <- decompose(nasa_ts)
plot(nasa_dts)


### Check for model accuracy

nasa_msts <- msts (df_monthly$value, seasonal.periods=c(12,120)) # #Yearly and Decade
nasa_msts_120_240 <- msts (df_monthly$value, seasonal.periods=c(12,120,240)) # # Yearly Decade and 2*Decade
nasa_msts_240 <- msts (df_monthly$value, seasonal.periods=c(12,240)) # # Yearly and 2*Decade


accuracy_df <- data.frame(model=c(""),runid=c(""),metric_mape=c(as.double("")),metric_mae=c(as.double("")))
# THE LOOP IS FOR ONLY MS_TS OBJECT MODELS
for (i in 1:10)
{ 
  nTest <- 12*i  
  nTrain <- length(nasa_msts)- nTest -1
  train <- window(nasa_msts, start=1, end=1+(nTrain)/(120))
  test <- window(nasa_msts, start=1+(nTrain+1)/(120), end=1+(nTrain+12)/(120))
  
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


## For sries 12, 240
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
  train <- ts(df_monthly$value,start = c(1880,1),end=c((2019-i),12),frequency = 12)
  test <- ts(df_monthly$value[length(train)+1:12],start = c((2019-i+1),1),frequency = 12)
  
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


#####
##### Rolling-horizon holdout: ARIMA on residuals
##### 

accuracy.arima=0 # we will check average 1-day-out accuracy for 7 days
for (i in 1:10)
{ 
  nTest <- 12*i  
  nTrain <- length(nasa_msts)- nTest -1
  train <- window(nasa_msts, start=1, end=1+(nTrain)/(120))
  test <- window(nasa_msts, start=1+(nTrain+1)/(120), end=1+(nTrain+12)/(120))
  
  cat("PREEEEE----------------------------------
      
      Data Partition",i,"
      
      Training Set includes",nTrain," time periods. Observations 1 to", nTrain, "
      Test Set includes 12 time periods. Observations", nTrain+1, "to", nTrain+12,"
      
      ")
  
  
  
  trainlm <- tslm(train ~ trend + season)
  trainlmf <- forecast(trainlm,h=12)
  
  residauto <- auto.arima(trainlm$residuals)
  residf <- forecast(residauto,h=12)
  
  y <- as.numeric(trainlmf$mean)
  x <- as.numeric(residf$mean)
  sp <- x+y
  
  cat("----------------------------------
      
      Data Partition",i,"
      
      Training Set includes",nTrain," time periods. Observations 1 to", nTrain, "
      Test Set includes 12 time periods. Observations", nTrain+1, "to", nTrain+12,"
      
      ")
  
  print(accuracy(sp,test))
  #  print(residauto)
  
  accuracy.arima<-rbind(accuracy.arima,accuracy(sp,test)[1,5])
  
  #print(sp$model)
}
accuracy.arima<-accuracy.arima[-1]

