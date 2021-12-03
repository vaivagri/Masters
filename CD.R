library(xts)
library(quantmod)
library(plyr)
library(tseries)
Sys.setlocale(locale = "English")

CD_full<-read.csv("CD.csv")
CD_full$Date<-as.Date(CD_full$Date)
CD_ts<-as.xts(CD_full$Close, order.by=as.Date(CD_full$Date), "%Y/%m/%d")
CD_df<-CD_full[, c(1,2)]
plot(CD_ts)

### Stationarity CD
# Pradiniai duomenys
print(paste('ADF=', adf.test(CD_ts)$p.value,',PP=', pp.test(CD_ts)$p.value, ',KPSS=', kpss.test(CD_ts)$p.value))
# Returns
CD_gr<-diff(CD_ts)[-1]
print(paste('ADF=', adf.test(CD_gr)$p.value, ',PP=', pp.test(CD_gr)$p.value, ',KPSS=', kpss.test(CD_gr)$p.value))
plot.xts(CD_gr)
# Log returns
CD_ret <- diff(log(CD_ts))[-1]    
print(paste('ADF=', adf.test(CD_ret)$p.value, ',PP=', pp.test(CD_ret)$p.value, ',KPSS=', kpss.test(CD_ret)$p.value))
plot.xts(CD_ret)

head(CD_ret)
summary(CD_ret)

### Stylized facts
# Volatility clusters
plot(CD_ret, ylab="", xlab="", main = "CD")
var1_ret100 <- tail(sort(abs(coredata(CD_ret))),200)[1]
idx <- which(coredata(CD_ret) <= var1_ret100) 
CD_retAbs100<-abs(CD_ret)
CD_retAbs100[idx]<-0
plot(CD_retAbs100)

# ACF
par(mfrow=c(1,2))
acf(coredata(CD_ret), main="ACF - CD") 
acf(abs(coredata(CD_ret)), main="ACF - abs. CD")
par(mfrow=c(1,1))

# Ljung-Box test
Box.test(CD_ret, type = "Ljung-Box", lag = 10) 
Box.test(abs(CD_ret), type = "Ljung-Box", lag = 10)

library(FinTS)
# LM test 
ArchTest(coredata(CD_ret), lags = 10) 

library(moments)
mean(CD_ret)
skewness(CD_ret)
kurtosis(CD_ret) 

# nonnormality and fat-tails (Leptokurtosis)
hist(CD_ret, probability = TRUE, main="Histogram of daily returns of CD")
curve(dnorm(x, mean=mean(CD_ret), sd=sd(CD_ret)), add=TRUE)
qqnorm(CD_ret)
qqline(CD_ret,col=2)
jarque.bera.test(as.vector(CD_ret)) 

############################################################################################################################################

library(rugarch)
library(tidyverse)
library(tidyr)
library(purrr)

#####
# Remiantis: PERLIN, M. S., M. MASTELLA, D. F. VANCIN, H. P. RAMOS. A GARCH tutorial with R. Revista de Administracao Contemporanea [interaktyvus]. 
# 2021, vol. 25(1), 1-16 [þiûrëta 2021-03-13]. e-ISSN 1982-7849. Prieiga per: https://doi.org/10.1590/1982-7849rac2021200088

max_lag_AR <- 1 
max_lag_MA <- 1 
max_lag_ARCH <- 2 
max_lag_GARCH <- 2
dist_to_use <- c('norm', 'std') 
models_to_estimate <- c('sGARCH', 'eGARCH')

library(dplyr)

out <- find_best_arch_model(x = CD_ret, 
                            type_models = models_to_estimate,
                            dist_to_use = dist_to_use,
                            max_lag_AR = max_lag_AR,
                            max_lag_MA = max_lag_MA,
                            max_lag_ARCH = max_lag_ARCH,
                            max_lag_GARCH = max_lag_GARCH)

tab_out <- out$tab_out

df_long <- tidyr::pivot_longer(data = tab_out %>%
                                 select(model_name,
                                        type_model,
                                        type_dist,
                                        AIC, BIC),  cols = c('AIC', 'BIC'))

models_names <- unique(df_long$model_name)
best_models <- c(tab_out$model_name[which.min(tab_out$AIC)],
                 tab_out$model_name[which.min(tab_out$BIC)])

# figure out where is the best model
df_long <- df_long %>%
  mutate(order_model = if_else(model_name %in% best_models, 'Best Model', 'Not Best Model') ) %>%
  na.omit()

# make table with best models
df_best_models <- df_long %>%
  group_by(name) %>%
  summarise(model_name = model_name[which.min(value)],
            value = value[which.min(value)],
            type_model = type_model[which.min(value)])

library(ggplot2)
# plot results
CD_p1 <- ggplot(df_long %>%
               arrange(type_model), 
             aes(x = reorder(model_name, 
                             order(type_model)),
                 y = value, 
                 shape = type_dist,
                 color = type_model)) + 
  geom_point(size = 3.5, alpha = 0.7) + 
  coord_flip() + 
  theme_bw(base_family = "TT Times New Roman") + 
  facet_wrap(~name, scales = 'free_x') + 
  geom_point(data = df_best_models, mapping = aes(x = reorder(model_name, 
                                                              order(type_model)),
                                                  y = value), 
             color = 'black', size = 4, shape = 8) +
  labs(title = 'GARCH modeliai - CD',
       x = '',
       y = 'AIC/BIC ivertis',
       shape = 'Skirstinio tipas',
       color = 'Modelio tipas') + 
  theme(legend.position = "bottom")

CD_p1

##### Best model
best_spec <- ugarchspec(variance.model=list(model="eGARCH", garchOrder=c(2,1)),mean.model=list(armaOrder=c(0,0),include.mean=TRUE),distribution.model="std")
fit_best <- ugarchfit(spec=best_spec,data=CD_ret)
fit_best

#####

par(mfrow=c(2,2))
plot(fit_best,which=6)
plot(fit_best,which=9)
plot(fit_best,which=10)
plot(fit_best,which=11)
par(mfrow=c(1,1))

####

# rolling  forecasting 
fit.best2 <- ugarchfit(spec=best_spec,data=CD_ret, out.sample=503) # (503 data for comparison)
(forc2 = ugarchforecast(fit.best2, n.ahead=1, n.roll=503)) 
plot(forc2, which=2) 
plot(forc2, which=4) 
