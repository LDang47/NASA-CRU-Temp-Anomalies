## Q1 and Q2 NASA
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


nasa_msts <- msts (df_monthly$value, seasonal.periods=c(12,120)) # #Yearly and Decade
nasa_msts_120_240 <- msts (df_monthly$value, seasonal.periods=c(12,120,240)) # # Yearly Decade and 2*Decade
nasa_msts_240 <- msts (df_monthly$value, seasonal.periods=c(12,240)) # # Yearly and 2*Decade

arima_fourier <- auto.arima(nasa_msts, seasonal=FALSE,
                            xreg=fourier(nasa_msts, K=c(5,5)))

plot(forecast(arima_fourier,xreg=fourier(nasa_msts, K=c(5,5), h=972)))

forecasted_values <- data.frame(forecast(arima_fourier,xreg=fourier(nasa_msts, K=c(5,5), h=972),level=c(80,95)))

forecasted_values$Date <- seq(as.Date('2020-01-01'),as.Date('2100-12-31'), "1 month")

write.csv(forecasted_values,"forecasted_nasa.csv")


checkresiduals(fit)
str(forecasted_values)

highchart(type = "stock") %>%
  hc_add_series(forecasted_values, "line", hcaes(Date,`Point.Forecast` ), name = "Forecasted", yAxis = 0) %>%
  hc_add_series(forecasted_values, "arearange", hcaes(Date, low = `Lo.95`, high = `Hi.95`), name = "Interval", yAxis = 0) %>% 
  hc_add_theme(hc_theme_ffx())%>%  hc_title(text = "<b>Global Mean Forecasting using ARIMA Fourier</b>",margin = 20, align =
                                              "center",style = list(color = "#4b4f57", useHTML = TRUE))


