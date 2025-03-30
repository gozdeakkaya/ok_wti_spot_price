### REQUIRED LIBRARIES
library(TSA)
library(forecast)
library(fpp2)
library(tseries)
library(dplyr)
library(lubridate)
library(gridExtra)
library(fUnitRoots)
library(caschrono)
library(ggplot2)
library(tidyquant)
library(tidyverse)
library(anomalize)
library(timetk)
library(knitr)
library(tibbletime)
library(lmtest)
library(uroot)
library(xts)
library(tsbox)
library(pdR)
library(tidyquant)
library(smooth)


### READING DATA
data <- read.table("monthly_data.txt", sep="," , header = F)
names(data) <- c("month","ok_wti_spot_price")
head(data)
dim(data)
data$month <- as.factor(data$month) #to use time_decompose() function, it is needed.

#Change Factor to Date format
data$month <- paste(data$month, "01", sep="-")
data$month <- as.Date(data$month,format="%Y-%m-%d")
df <- as_tibble(data)
class(df)

#Anomaly Detection
df_anomalized <- df %>%
  time_decompose(ok_wti_spot_price, merge = TRUE) %>%
  anomalize(remainder) %>%
  time_recompose()
df_anomalized %>% glimpse()

sum(df_anomalized$anomaly == "Yes") #There is 30 anomalies in the data.

#Anomaly Points Visualization
anomaly_plot <- df_anomalized %>%
  plot_anomaly_decomposition() +
  ggtitle("Freq/Trend = 'auto'")
anomaly_plot

#Cleaning the anomalies from the data
data_tibble_clean <- df %>%
  time_decompose(ok_wti_spot_price, merge = TRUE) %>%
  anomalize(remainder) %>%
  clean_anomalies() %>%
  time_recompose()
data_tibble_clean %>% glimpse()
sum(data_tibble_clean$anomaly == "Yes")

#Showing change in observed versus observed_clean
data_tibble_clean %>% 
  filter(anomaly == "Yes") %>%
  dplyr::select(month, anomaly, observed, observed_cleaned) %>%
  kable()
class(data_tibble_clean)

ts <- as.data.frame(data_tibble_clean) #changing tibble to data frame

#changing data.frame to ts class
x <- data_tibble_clean %>%
  dplyr::select(month, observed_cleaned) %>%
  as.ts()

data_ts <- ts(x[,2],start=c(1986,1),frequency=12)
data_ts

class(data_ts)

### DATA DESCRIPTION
summary(data_ts)
autoplot(data_ts)

### SPLITTING THE DATA INTO TRAIN AND TEST SET
#view(data)
train <- window(data_ts, end = c(2020,3)) #411 obs
test <- window(data_ts, start = c(2020,4)) #12 obs

data_ts_trans <- log(train)
boxcox(lm(train ~ 1))

plot(data_ts_trans, ylab='ok_wti_spot_price', xlab='month', type='o', main="Time Series Plot of WTI Crude Oil Spot Price")
p1 <- ggAcf(data_ts_trans,main="ACF of WTI Crude Oil Spot Price")
p2 <- ggPacf(data_ts_trans,main="PACF of WTI Crude Oil Spot Price")
grid.arrange(p1,p2,nrow=1)

### FORMAL TESTS

#Checking stationary, trend and seasonality by conducting formal tests.
#KPSS Test
kpss.test(data_ts_trans,null=c("Level")) 
kpss.test(data_ts_trans,null=c("Trend"))

#ADF Test
mean(data_ts_trans) #There is a constant.
adfTest(data_ts_trans, lags=2, type="c")  #ADF test shows the process is non stationarity.
adfTest(data_ts_trans, lags=2, type="ct") #ADF test shows that there is a stochastic trend.

#HEGY Test
HEGY.test(wts=data_ts_trans, itsd=c(1,0,0), regvar=0, selectlags=list(mode="signf", Pmax=NULL))$stats
ndiffs(data_ts_trans)
nsdiffs(diff(data_ts_trans))

#Making the process stationary by taking one regular differencing
data_ts_trans_diff <- diff(data_ts_trans)
autoplot(data_ts_trans_diff, main="Time Series Plot of Differenced Series")
#The time series plot shows that the process became stationary after differencing.

p1 <- ggAcf(data_ts_trans_diff,main="ACF Plot of  Stationary Data Set")
p2 <- ggPacf(data_ts_trans_diff,main="PACF Plot of  Stationary Data Set")
grid.arrange(p1,p2,nrow=1)

#I have checked the stationarity of the series one more time by conducting formal tests.
kpss.test(data_ts_trans_diff, null=c("Level"))
mean(data_ts_trans_diff)
adfTest(data_ts_trans_diff, lags=1, type="nc")
ndiffs(data_ts_trans_diff)
nsdiffs(data_ts_trans_diff)
#After taking one regular differencing, the time series became stationary according to the formal tests.


### MODEL SELECTION
#Identifying a proper ARIMA model.
#Using ESACF method ARIMA(0,1,0), ARIMA(1,1,1)
eacf(data_ts_trans_diff) 
fit1 <- Arima(data_ts_trans, order=c(1,1,1),seasonal = c(0,0,0)) #NOT SIGNIFICANT
fit1

#Using armaselect function
armaselect(data_ts_trans_diff)
fit2 <- Arima(data_ts_trans,order=c(0,1,1),seasonal = c(0,0,0))  #SIGNIFICANT
fit2  #-845.86

fit3 <- Arima(data_ts_trans, order=c(0,1,2),seasonal = c(0,0,0)) #NOT SIGNIFICANT
fit3  #-845.23 

fit4 <- Arima(data_ts_trans,order=c(1,1,2),seasonal = c(0,0,0)) #SIGNIFICANT
fit4  #-848.84

fit5 <- Arima(data_ts_trans, order=c(2,1,1),seasonal = c(0,0,0)) #SIGNIFICANT
fit5  #-849.91

fit6 <- Arima(data_ts_trans, order=c(1,1,0),seasonal = c(0,0,0)) #SIGNIFICANT
fit6  #-845.03

fit7 <- Arima(data_ts_trans, order=c(0,1,3),seasonal = c(0,0,0)) #SIGNIFICANT
fit7  #-847.44

#I have chosen the model which has minimum AIC value. ARIMA(2,1,1) fit5

### DIAGNOSTIC CHECKING

#Normality of Residuals
r <- resid(fit5)
autoplot(r)+geom_line(y=0)+theme_minimal()+ggtitle("Plot of The Residuals")
ggplot(r, aes(sample = r)) +stat_qq()+geom_qq_line()+ggtitle("QQ Plot of the Residuals")+theme_minimal()
ggplot(r,aes(x=r))+geom_histogram(bins=20)+geom_density()+ggtitle("Histogram of Residuals")+theme_minimal()

jarque.bera.test(r)
shapiro.test(r) 


#Detection of Serial Autoorrelation
ggAcf(as.vector(r), main="ACF of the Residuals", lag = 48)+theme_minimal() 
#In order to say that we do not have any correlation problem, all spikes should be in the White Noise band.
#There are some spikes that outside in the WN band.

m <- lm(r ~ 1+zlag(r))
bgtest(m)
Box.test(r, lag=2, type = c("Ljung-Box")) #Since p value is greater than alpha, we have 95% confident that the residuals of the model are uncorrelated, according to results of Box-Ljung Test.
Box.test(r, lag=2, type = c("Box-Pierce")) # Since p value is greater than alpha, we have 95% confident that the residuals of the model are uncorrelated, according to results of Box-Pierce Test.

#Heteroscedasticity
rr <- r^2
g1 <- ggAcf(as.vector(rr))+theme_minimal()+ggtitle("ACF of Squared Residuals")
g2 <- ggPacf(as.vector(rr))+theme_minimal()+ggtitle("PACF of Squared Residuals")
grid.arrange(g1,g2,ncol=2) 
#It is clearly seen that the all squared residuals are not in the 95% White Noise Band. Therefore, it is said that the errors are not homoscedastic.
#Both plots shows that some spikes are out of the white noise bands that is an indication of heteroscedasticity problem.

m <- lm(r ~ data_ts_trans+zlag(data_ts_trans)+zlag(data_ts_trans,2)) #White test 
bptest(m) # Since p value is greater than alpha, we fail reject Ho. Therefore, we can say that we have enough evidence to claim that there is no heteroscedasticity problem, according to results of Breusch-Pagan test.

### FORECASTING

#Performing Minimum MSE Forecast for the stochastic models 
#ARIMA
forecast1 <- forecast(fit5, h=12)
autoplot(forecast1)+theme_minimal()+ggtitle("Forecast of ARIMA")
ff <- exp(forecast1$mean)
shapiro.test(forecast1$residuals) #not normal

accuracy(ff,test)

#ETS
ets_forecast <- ets(train, model="MNN") #The algorithm automatically decides type of the components.
ets_forecast
fr1 <- forecast(ets_forecast,h=12) 
fr1
autoplot(fr1)+autolayer(fitted(fr1),series="fitted")+theme_minimal()
shapiro.test(fr1$residuals) #not normal



#NEURAL networks
nnmodel <- nnetar(train)
nnmodel
autoplot(train)+autolayer(fitted(nnmodel))+theme_minimal()+ggtitle("Fitted Values of Neural Network Model")
nnforecast <- forecast(nnmodel,h=12,PI=T)
nnforecast
autoplot(nnforecast)+theme_minimal()

shapiro.test(nnmodel$residuals) #not normal

accuracy(nnforecast,test)
accuracy(fr1,test)
accuracy(ff,test)




