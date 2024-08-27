# ---------------------------------------------------------------------------- #
# Estimates of Output Gap for Colombia
# Methodologies included in the code are: 
# - HP Filter 
# - Never use HP Filter (Hamilton)
# - Kalman Filter
# - Production Function Approach (includes NAIRU and NAICU estimates)
# 
#
# - Author: Juan Camilo Meneses Cortés - Universidad de los Andes
# - Version: 1.0 
# - Last update: Jan-2021
# ---------------------------------------------------------------------------- #

#Libraries 
library(readxl)
library(tidyr)

library(seasonal)
library(splines2)
library(splines)
library(mgcv)
library(rms)

library(KFAS)
library(dlm)
library(dse)

library(mFilter)
library(neverhpfilter)
library(forecast)

rm(list = ls())
#Working Directory
setwd("C:/Users/juanm/OneDrive - Universidad de los Andes/Juan Camilo Meneses/Personal/Tesis PEG/Datos")

#Read Data: 
data = read_excel("Datos.xlsx", sheet = "Data", col_names = T, skip = 0)
names(data)

   # ts data type 
   data_ts = ts(data[,-1], start = c(1994,1), frequency = 4)

   # xts data type
   data_xts = xts(data[,-1], as.Date(data$Fecha))

#Auxiliary functions
   reverse_ts <- function(y)
   {
     ts(rev(y), start=tsp(y)[1L], frequency=frequency(y))
   }
   
   reverse_forecast <- function(object)
   {
     h <- length(object[["mean"]])
     f <- frequency(object[["mean"]])
     object[["x"]] <- reverse_ts(object[["x"]])
     object[["mean"]] <- ts(rev(object[["mean"]]),
                            end=tsp(object[["x"]])[1L]-1/f, frequency=f)
     object[["lower"]] <- object[["lower"]][h:1L,]
     object[["upper"]] <- object[["upper"]][h:1L,]
     return(object)
   }

# ---------------------- Adjustments (SA - Filters) --------------------------- 
# --------------------- For additional Series -------------------------------- #
   
#Seasonal Adjustment of Inflation: 
 cpi_core_seas = ts(data$inflation_a_cpi_core, start = c(1994,1), frequency = 4) 
 cpi_core_seasx13 = seas(cpi_core_seas)
 plot(cpi_core_seasx13)
 
 cpi_core_seas_15 = ts(data$inflation_core_15, start = c(1994,1), frequency = 4) 
 cpi_core_seasx13_15 = seas(cpi_core_seas_15)
 plot(cpi_core_seasx13_15)
 
 cpi_core_sa = as.data.frame(final(cpi_core_seasx13))
 write.csv(cpi_core_sa, "cpi_core_sa.csv")

#Seasonal Adjustment of historical (:2000) unemployment: 
 historical_unemployment = read_excel("Datos.xlsx", sheet = "NAIRU_Unemployment", 
                           col_names = T, skip = 0, range = "A2:D64")
 historical_unemployment = ts(historical_unemployment[,-1], start = c(1984,4),
                              frequency = 4)
 historical_unemployment_seas = seas(historical_unemployment[,1])
 plot(historical_unemployment_seas)
 historical_unemployment_sa = as.data.frame(final(historical_unemployment_seas))
 
 historical_employment_seas = seas(historical_unemployment[,2])
 plot(historical_employment_seas)
 historical_employment_sa = as.data.frame(final(historical_employment_seas))
 
 historical_ls_seas = seas(historical_unemployment[,3])
 plot(historical_ls_seas)
 historical_ls_seas = as.data.frame(final(historical_ls_seas))
 
 
 write.csv(cbind(historical_unemployment_sa, historical_employment_sa,
                 historical_ls_seas), "historical_unemployment_sa.csv")

#Seasonal Adjustment of Installed Capacity Utilization:
 plot(data$capacity_u, type = "l")
 icu = ts(data$capacity_u, start = c(1994,4), frequency = 4)
 icu_stl = stl(icu, s.window = "periodic", t.window = 13, robust = TRUE)
 plot(icu_stl)
 icu_trend_cycle = icu_stl$time.series[,2]
 plot(icu_trend_cycle)
  #H-P Filter Adjustment over ICU trend-cycle series: 
  icu_hp = hpfilter(icu_trend_cycle, type = "lambda", freq = 1600)$trend
  plot(icu_hp, type = 'l', col = "blue")
  lines(icu_trend_cycle)
  write.csv(cbind(icu_trend_cycle, icu_hp), "icu_decomposed.csv")
  
# Filter for the real exchange rate GAP  
  rer_stl = stl(ts(data$real_exchange_rate, start = c(1994,1), frequency = 4), t.window=13, s.window="periodic", robust=TRUE)
  plot(rer_stl$time.series[,2])
  rer_trend_cycle = rer_stl$time.series[,2]
  
  rer_hp = hpfilter(rer_trend_cycle, freq = 1600, type = c("lambda"))
  rer_trend = rer_hp$trend
  plot(data_xts$real_exchange_rate)
  plot(rer_trend, col = 'blue')
  lines(rer_trend_cycle)
  lines(ts(data$real_exchange_rate, start = c(1994,1), frequency = 4))
  write.csv(cbind(rer_trend_cycle, rer_trend), "rer_trend.csv")  
  
# Filter for the real exchange rate depreciation  
  
  rer_dep_hp = hpfilter(na.omit(ts(data_xts$rer_depreciation, start = c(1994,1), frequency = 4)), freq = 1600, type = c("lambda"))
  rer_dep_hp = rer_dep_hp$trend
  plot(ts(data_xts$rer_depreciation, start = c(1994,1), frequency = 4))
  lines(rer_dep_hp, col = 'blue')
  write.csv(rer_dep_hp, "rer_dep_trend.csv")    
  
# Filter for EMBI col
  plot(data$EMBI_col, type = "l")
  EMBI_stl = stl(na.omit(ts(data$EMBI_col, start = c(1994,1), frequency = 4)), t.window=13, s.window="periodic", robust=TRUE)
  EMBI_trend_cycle = EMBI_stl$time.series[,2]
  plot(EMBI_trend_cycle)
  
  EMBI_hp = hpfilter(EMBI_trend_cycle, freq = 1600, type = c("lambda"))
  EMBI_hp = EMBI_hp$trend
  plot(EMBI_trend_cycle)
  lines(EMBI_hp, col = 'blue')
  write.csv(cbind(EMBI_trend_cycle,EMBI_hp), "EMBI_decomposed.csv")    

# ----------------- Univariate Filters for Outputgap ---------------------------   
# ---------------------------- HP Filter ------------------------------------- #
# Lambda parameter is set to 1600 for quarterly data: 
 
 #simple H_P Filter   
 outputgap_hp_1600 = as.ts(hpfilter(data_ts[,2], type = "lambda", freq = 1600)$trend)
 plot(outputgap_hp_1600)
 outputgap_hp_100 = as.ts(hpfilter(data_ts[,2], type = "lambda", freq = 100)$trend)
 plot(outputgap_hp_100)
 write.csv(cbind(outputgap_hp_100, outputgap_hp_1600), "outputgap_hp.csv")

 #Extended HP-Filter: 
  #Forecast and Backcast series of GDP annual growth (Julio 2011): 
  gdp_annual = diff(data_ts[1:eval(nrow(data_ts)-2),2],4) #1995Q1 - 2019Q4
  gdp_annual_reverse = rev(gdp_annual)
  
  plot(gdp_annual, type = "l")
  model_gdp_annual = arima(gdp_annual, order = c(2,0,4), method = "ML")
    #Diagnosis Test
    Box.test(residuals(model_gdp_annual), lag = nrow(data_ts)/4)
    jarque.bera.test(residuals(model_gdp_annual))
    checkresiduals(model_gdp_annual)
    acf2(residuals(model_gdp_annual))
    plot(forecast(model_gdp_annual, h = 20))
    forecast_gdp = forecast(model_gdp_annual, h = 20)$mean
    
  model_gdp_annual_reverse = arima(gdp_annual_reverse, order = c(2,0,4), method = "ML")
    #Diagnosis Test
    Box.test(residuals(model_gdp_annual_reverse), lag = nrow(data_ts)/4)
    jarque.bera.test(residuals(model_gdp_annual_reverse))
    checkresiduals(model_gdp_annual_reverse)
    acf2(residuals(model_gdp_annual_reverse))  
    plot(reverse_forecast(forecast(model_gdp_annual_reverse, h = 20)))  
    backcast_gdp = reverse_forecast(forecast(model_gdp_annual_reverse, h = 20))$mean
    
    
 gdp_annual = as.data.frame(gdp_annual)     
 colnames(gdp_annual) = c("x")
 gdp_annual_extended = c(t(as.vector(backcast_gdp)), t(as.vector(gdp_annual)),as.vector(forecast_gdp))
 gdp_annual_extended = ts(gdp_annual_extended, start = c(1989,1), frequency = 4)
 plot(gdp_annual_extended)
 write.csv(gdp_annual_extended, "GDP_Annual_Var_Extended.csv")
 
 gdp_extended = ts(read_excel("Datos.xlsx", sheet = "Output Gap", col_names = T, skip = 0, range = "C1:C141"), start = c(1990,1), frequency = 4)
 gdp_extended = ts(as.vector(log(gdp_extended)), start = c(1990,1), frequency = 4)
 stl_lgdp = stl(gdp_extended, s.window = 13, robust = TRUE)
 plot(stl_lgdp)
 lgdp_trend_cycle = stl_lgdp$time.series[,2]
 write.csv(lgdp_trend_cycle, "lgdp_trend_cycle.csv")
 
 #Hp Filter over Trend - Cycle Component: 
 outputgap_hp__extended_1600 = as.ts(hpfilter(lgdp_trend_cycle, type = "lambda", freq = 1600)$trend)
 plot(outputgap_hp__extended_1600)
 write.csv(outputgap_hp__extended_1600, "outputgap_hp_extended.csv")
 
# ------------------------. Never  HP Filter --------------------------------- #
# Look-ahead period for quarterly data = 8. Lags (seasonality) default = 4:     
 gdp_extended_xts = read_excel("Datos.xlsx", sheet = "Output Gap", col_names = T, skip = 0, range = "A1:F141")
 gdp_extended_xts = xts(gdp_extended_xts[,-1], as.Date(gdp_extended_xts$Date))
 lgdp_trend_cycle_xts = gdp_extended_xts$LGDP_T_C["1990-03-01/2024-12-01"]
 outputgap_nhp_1qr = yth_filter(lgdp_trend_cycle_xts, h = 6, p = 4, output = "trend")
 plot(outputgap_nhp_1qr)
 plot(lgdp_trend_cycle_xts - outputgap_nhp_1qr)
 
 lgdp_trend_nhp = NULL
 for(i in 2:14){
   aux = yth_filter(lgdp_trend_cycle_xts, h = i, p = 4, output = "trend")
   lgdp_trend_nhp = cbind(lgdp_trend_nhp, aux)
 }
 lgdp_trend_nhp = na.omit(lgdp_trend_nhp)
 lgdp_trend_nhp_final = rowMeans(lgdp_trend_nhp)
 lgdp_trend_nhp_final = hpfilter(as.ts(lgdp_trend_nhp_final), type = "lambda", freq = 1600)$trend
 plot(lgdp_trend_nhp_final)
 write.csv(lgdp_trend_nhp_final, "outputgap_nhp.csv", row.names = T)

# ----------------- Structural Approaches for Outputgap ------------------------ 
# ------------------------. Production Function ------------------------------ # 
# Cobb - Douglas Production Function
  #NAIRU estimates:
  #Based on Trend-cycle component: 
  unemployment_stl = stl(ts(data$Unemployment_u, start = c(1994,1), frequency = 4),t.window=13, s.window="periodic", robust=TRUE) 
  plot(unemployment_stl)
  unemployment_trend_cycle = xts(unemployment_stl$time.series[,2], as.Date(data$Fecha)) 
  
  data_xts = cbind(data_xts, unemployment_trend_cycle)
  plot(data_xts$unemployment_trend_cycle)
  
   data_nairu = cbind.xts(diff(data_xts$inflation_cpi_core),
            lag(diff(data_xts$inflation_cpi_core)), lag(diff(data_xts$inflation_cpi_core,2))
         ,lag(diff(data_xts$inflation_cpi_core,3)), lag(diff(data_xts$inflation_cpi_core,4))
         ,data_xts$unemployment_trend_cycle, lag(data_xts$unemployment_trend_cycle,1), 
         lag(data_xts$unemployment_trend_cycle,2), lag(data_xts$unemployment_trend_cycle,3), 
         lag(data_xts$unemployment_trend_cycle,4), lag(data_xts$inflation_food,1), lag(data_xts$import_inflation_ppi,1))
    data_nairu = na.omit(data_nairu)
    names(data_nairu)
    spline_data = seq(1,nrow(data_nairu),1)
    spline_data_2 = spline_data^2
    spline_data_3 = spline_data^3
  
  data_nairu = cbind.xts(data_nairu, spline_data, spline_data_2,spline_data_3)
  
  colnames(data_nairu) = c("inflation", "inflation.l1", "inflation.l2",
                           "inflation.l3",  "inflation.l4","unemployment", "unemployment.l1"
                           , "unemployment.l2", "unemployment.l3", "unemployment.l4"
                           , "inflation_food", "inflation_import", "sp1", "sp2", "sp3")
  
  #Lm estimation of Phillips Curve and spline variables:
  phillips_curve = lm(inflation ~ 1 + sp1 + sp2 + sp3 + inflation.l1+inflation.l2 + inflation.l3 + inflation.l4
                      +unemployment+unemployment.l1 + unemployment.l2+ unemployment.l3 + unemployment.l4+ 
                     inflation_food+inflation_import,data = data_nairu)
  summary(phillips_curve)
  phi = c(coefficients(phillips_curve))[c(1:4)]
  b1 = sum(coefficients(phillips_curve)[c(9,10,11,12, 13)])
  nairu_spline_0 = NULL
  for(i in 1: nrow(data_nairu)){
     nairu_spline_0[i] = -(1*phi[1] + phi[2]*data_nairu$sp1[i,]+phi[3]*data_nairu$sp2[i,]+phi[4]*data_nairu$sp3[i,])/b1
  }
  plot(nairu_spline_0)
  
  #Lm estimation of Cubic Spline with 3 knots:
  phillips_curve_m = lm(inflation ~ bs(sp1, knots = c(15,60,85))+1+inflation.l1+inflation.l2 + inflation.l3 + inflation.l4
                         +unemployment+unemployment.l1 + unemployment.l2+unemployment.l3+unemployment.l4+ 
                        inflation_food+inflation_import,data = data_nairu)
  summary(phillips_curve_m)
  phi = c(coefficients(phillips_curve_m))[c(1:7)]
  b1 = sum(coefficients(phillips_curve_m)[c(12,13,14,15,16)])
  spline_terms = phillips_curve_m$model$`bs(sp1, knots = c(15, 60, 85))`
  nairu_spline_1 = NULL
  for(i in 1: nrow(spline_terms)){
     nairu_spline_1[i] = -(1*phi[1] + phi[2]*spline_terms[i,1]+phi[3]*spline_terms[i,2]+phi[4]*spline_terms[i,3]+phi[5]*spline_terms[i,4]+phi[6]*spline_terms[i,5]+phi[7]*spline_terms[i,6])/b1
  }
  plot(nairu_spline_1)
  lines(phillips_curve_m$model$unemployment)
  
  #constant NAIRU:
  summary(lm(inflation ~ 1+ inflation.l1+inflation.l2 + inflation.l3 + inflation.l4
     +unemployment+unemployment.l1 + unemployment.l2+ unemployment.l3 + unemployment.l4+ 
        inflation_food+inflation_import,data = data_nairu))
  lm_constant_NIARU = lm(inflation ~ 1+ inflation.l1+inflation.l2 + inflation.l3 + inflation.l4
     +unemployment+unemployment.l1 + unemployment.l2+ unemployment.l3 + unemployment.l4+ 
       inflation_food+inflation_import,data = data_nairu)
  coef_constant_nairu = coef(lm_constant_NIARU)
  nairu_constant = -(coef_constant_nairu[1])/(sum(coef_constant_nairu[6:10])) 
 
  #Restricted Cubic Spline Just With Unemployment: 
  unemployment_curve = lm(unemployment ~ 1 + bs(sp1,knots = c(40,80)),data = data_nairu)
  summary(unemployment_curve)
  nairu_spline_u = predict(unemployment_curve, data_nairu)
  plot(nairu_spline_u)
  lines(unemployment_curve$model$unemployment)
  
  nairu_hp = hpfilter(as.ts(data_nairu$unemployment), type = "lambda", freq = 1600)$trend
  plot(nairu_hp)
  lines(as.ts(data_nairu$unemployment))
  write.csv(cbind(nairu_hp, nairu_spline_u, nairu_spline_0, nairu_spline_1, nairu_constant, data_nairu$unemployment), "NAIRU.csv")

 # ------------------------- Kalman Filter for NAICU  ------------------------ # 
  # DLM Package:
  #Define function for MLE estimation: 
  my_dlmfc <- function(par = c(phi_1, phi_2, phi_3, phi_4, phi_5,phi_6, phi_7,gamma)){
     phi_1 = par[1]
     phi_2 = par[2]
     phi_3 = par[3]
     phi_4 = par[4]
     phi_5 = par[5]
     phi_6 = par[6]
     phi_7 = par[7]
     gamma = par[8]
     
     #Transition matrix G
     GG =  matrix(c(phi_1,0,gamma,0,
                    0,1,0,0,
                    0,0,phi_2, phi_3,
                    0,0,1,0
                    ), nrow = 4, byrow = T)
     
     #Measurement equation matrix
     FF =   matrix(c(0,1,1,0,
                    1,0,0,0), nrow = 2, byrow = T)
     
     #Measurement equation var-cov matrix: 
     V = diag(0, 2,2)
     
     #State equation var-cov matrix: 
     W  = matrix(c(phi_4, 0,0,0,
                   0,phi_5,0,0,
                   0,0,phi_6,0,
                   0,0,0,phi_7), nrow = 4, byrow = T)
     
     m0 = matrix(c(0,0,0,0), nrow = 4)
     C0 = matrix(rep(0,16), nrow = 4)
     my_dlm1 = dlm(FF=FF,V=V,GG=GG,W=t(W)*W,V=V,m0=m0,C0=C0)
  }
  #Data for estimation (capacity utilization and inflation): 
  data_naicu = cbind(data$capacity_u_tc,data$inflation_cpi_core)
  
  #f (!(all.equal(W_prueba, t(W_prueba)) && all(eigen(W_prueba)$values >= 
                                                 0))) 
  #stop("W is not a valid variance matrix")
  
  #parm = c(0.86, 0.97, -0.001, 1, 1, 0.07,1, -9)
  my_dlm__naiuc_mle = dlmMLE(data_naicu, parm = c(0.8, 0.9, -0.001, 1, 1,0.07,1,-10),lower=c(0,0,0,-Inf,-Inf -Inf,-Inf,0), my_dlmfc)
  my_dlm__naiuc_mle$convergence
  print(my_dlm__naiuc_mle)
  model.fit = my_dlmfc(my_dlm__naiuc_mle$par)
  print(model.fit)
  my_dlm_smoothed <- dlmSmooth(data_naicu, model.fit) 
  my_dlm_smoothed <- my_dlm_smoothed$s
  plot(my_dlm_smoothed)
  #NAICU 
  naicu_dlm_1 = my_dlm_smoothed[,2]
  plot(naicu_dlm_1, type = "l")
  lines(data_naicu[,1])
  
  model.build <- function(parm){
     return(
        dlmModPoly(order = 4, dV = exp(parm[1]), dW = exp(parm[2:5]))
     +dlmModSeas(4, dV=parm[6])
     )
  }
  
  model.mle <- dlmMLE(data_naicu[,1], parm = c(rep(0,6)), build=model.build, hessian = TRUE)
  if(model.mle$convergence==0) print("converged") else print("did not converge")
  model.mle$par 
  model.fit <- model.build(model.mle$par)
  model.filtered <- dlmFilter(data_naicu[,1], model.fit)
  model.smoothed <- dlmSmooth(model.filtered)$s
  plot(model.smoothed[,1])
  #naicu = apply(model.smoothed[-1,1:5],1, sum)
  naicu = model.smoothed[,1]
  windows()
  plot(naicu, type = 'l', col = 'seagreen')
  lines(data_naicu[,1])
  
  n <- 1*4
  model.forecast <- dlmForecast(model.filtered, nAhead=n)
  
  x <- index(data_naicu)
  xf <- seq(max(x), max(x)+n/4, 1/4)
  aa <- model.forecast$a[,-1]*(-1)
  aa <- cbind(model.forecast$a[,1], aa)
  a <- drop(model.forecast$a%*%t(FF(model.fit)))
  a <- c(tail(data_naicu,1), a)
  df <- rbind(
     data.frame(x=x, y=as.numeric(data_naicu), series="original"),
     data.frame(x=x, y=model.filtered$m[-1,1], series="filtered"),
     data.frame(x=x, y=model.smoothed[-1,1], series="smoothed"),
     data.frame(x=xf, y=a, series="forecast")
  )
  g.dlm <- ggplot(subset(df, x>10), aes(x=x, y=y, colour=series)) + geom_line()
  windows()
  g.dlm
  write.csv(naicu, 'naicu_dic.csv')
  
  
  # ----------- Kalman Filter with dse package --------------------- #
  #Define matrices: 
  #Transition matrix F (state equation): 
  
  f = matrix(c(0.86,0,0.007,0,
               0,1,0,0,
               0,0,0.5, -0.1,
               0,0,1,0
  ), nrow = 4, byrow = T)
  
  # System noise distribution matrix Q: 
  q = diag(1.1, 4,4) #exp(1) is selected firs
  #q[1,1] <- 2.567482e-05
  #q[2,2] <- 0.85*8.696542e-05
  #q[3,3] <- 8.696542e-05
  q[4,4] <- 0
  #Input matrix G (state equation, exogenous): 
  g = matrix(c(0.05,
               0,
               0,
               0), nrow = 4, byrow = T)  
  #Output matrix H (measurement equation): 
  h = matrix(c(0,1,1,0,
               1,0,0,0), nrow = 2, byrow = T)
  #Output measurement noise distribution matrix R:
  r = diag(0, 2,2)
  
  
  #Construct state-space model (no innovations form):
  model_k1_icu= SS(F.= f, G = g, H = h, Q = q, R = r, z0 = c(0.7,0.03))
  data_dse_icu = TSdata(input = data$import_inflation_ppi[5:eval(nrow(data)-2)], output = cbind(data$capacity_u_tc, data$inflation_cpi_core)[7:nrow(data_ts),])
  windows()
  tfplot(data_dse_icu)
  #Estimation of model by maximum likelihood:
  model_k1_ml_icu = estMaxLik(model_k1_icu, data_dse_icu, algorithm = "optim",
                              algorithm.args=list(method="L-BFGS-B", upper=c(1,Inf,0.98,0.98,0.98,Inf, Inf, Inf), lower=c(0,0,-Inf, -Inf,0,0,0,0), hessian=TRUE,
                                                  control=list(maxit=10000))) 
  summary(model_k1_ml_icu)
  print(model_k1_ml_icu)
  
  #Smoother for obtaining the state filtered values: 
  smoother_model1_icu = smoother(model_k1_ml_icu)
  #Potential output and output GAP:
  naicu_k1 = smoother_model1_icu$smooth$state[,2]
  windows()
  plot(naicu_k1, type='l',col = 'seagreen')
  lines(data_dse_icu$output[,1])
  plot(smoother_model1_icu$smooth$state[,2], type = 'l', col = 'blue')
  lines(data_dse_icu$output[,1], type = 'l')
  write.csv(as.data.frame(naicu_k1), "naicu_kf1.csv")
  
  
# OLS regression for Solow Residual: 
  names(data)
  #Decomposition of series: 
  employment_stl = stl(ts(data$employment, start = c(1994,1), frequency = 4),  s.window = "periodic", t.window = 13, robust = TRUE)
  employment_trend_cycle = employment_stl$time.series[,2]
  plot(employment_trend_cycle)
  
  labor_stl = stl(ts(data$labor_supply, start = c(1994,1), frequency = 4),  s.window = "periodic", t.window = 13, robust = TRUE)
  labor_trend_cycle = labor_stl$time.series[,2]
  plot(labor_trend_cycle)
  
  write.csv(cbind(employment_trend_cycle, labor_trend_cycle), "employ_labors_tc.csv")
  capital_stock_stl = stl(na.omit(ts(data$capital_stock, start = c(1994,1), frequency = 4)),  s.window = "periodic", t.window = 13, robust = TRUE)
  capital_stock_trend_cycle = capital_stock_stl$time.series[,2]
  plot(capital_stock_stl)
  
  write.csv(capital_stock_trend_cycle, "capital_stock_tc.csv")
  
  #Data for CB function
  data_cb = cbind.xts(data_xts$lgdp_tc, data_xts$labor_supply_tc,data_xts$employment_tc,
                      data_xts$nairu,data_xts$capital_stock_tc,
                      data$capacity_u_tc, data_xts$naicu_promedio)
  capital_icu = lag(data_cb$capital_stock)*data_cb$data.capacity_u
  colnames(capital_icu) = c("capital_stock_icu")
  data_cb = na.omit(cbind.xts(data_cb, capital_icu))
  plot(data_cb$lgdp_tc)
  
  cb_ols = lm(lgdp_tc ~ log(capital_stock_icu) + log(employment_tc)-1, data = data_cb)
  summary(cb_ols)
  plot(residuals(cb_ols), type = 'l')
  cb_ols_res = lm(lgdp_tc ~ I(log(employment_tc)-log(capital_stock_icu)) + offset(log(capital_stock_icu))-1, data = data_cb)
  summary(cb_ols_res)
  library(limSolve)
  lsei(cbind.xts(log(data_cb$capital_stock_icu),log(data_cb$employment_tc)), data_cb$lgdp_tc, c(1, 1), 1)
  
  #Solow residual 
  ptf = data_cb$lgdp_tc - 0.766595*log(data_cb$employment_tc) - (1- 0.766595)*log(data_cb$capital_stock_icu)
  colnames(ptf) = c("total_factor_productivity")
  plot(ptf)
  
  ptf_ts = ts(ptf, start = c(1998,2), frequency = 4)
  ptf_trend = hpfilter(ptf_ts, type = "lambda", freq = 1600)$trend
  plot(ptf_trend, type = 'l')
  lines(ptf_ts)
  ptf_trend = as.xts(ptf_trend)
  write.csv(cbind(ptf_ts, ptf_trend), "tfp.csv")
  
  #Potential GDP from Cobb-Douglas Production Function: 
  y_pot_log = as.data.frame(ptf_trend) + (1-0.76603)*(lag(log(data_cb$capital_stock))+log(data_cb$naicu)) + 0.76603*(log(data_cb$labor_supply)+ log(1-data_cb$nairu))
  

# ---------------------------------------------------------------------------- #    
# ------------------ Kalman Filter with the dse package ---------------------- #
# ------------------------ Estimates for Output GAP -------------------------- #  

data_kalman = cbind(data_xts$inflation_core_15, data_xts$import_inflation_ppi,
                    data_xts$rer_depreciation, data_xts$lgdp_tc, data_xts$repo_real)["1999-03-01/2020-01-01"] 
prono_kalman = NULL
for(i in 1:ncol(data_kalman)){
  arima_aux = auto.arima(data_kalman[,1])
  forecast_aux = forecast(arima_aux, h = 20)$mean
  prono_kalman = cbind(prono_kalman, forecast_aux)
}
data_kalman_extended = rbind(as.data.frame(data_kalman, prono_kalman))
  
data_kalman_ts = ts(data_kalman, start = c(1999,1), frequency = 4)  

#Define matrices: 
 #Transition matrix F (state equation): 
  
  f = matrix(c(1,1,0,0,0,0,0,
               0, 0.8, 0.2,0,0,0,0,
               0,0,1,0,0,0,0,
               0,0,0,0.8,0,-0.99,0,
               0,0,0,0,1,0,0,
               0,0,0,0,0,0,0,
               0,0,0,0.4,0,0,0.8), nrow = 7, byrow = T)
 # System noise distribution matrix Q: 
 q = diag(0.5, 7,7)
 q[3,3] <- 0
 #Input matrix G (state equation, exogenous): 
 g = matrix(c(0,0,
              0,0,
              0,0,
              0.5,0,
              0,0,
              0,0,
              0,0.5), nrow = 7, byrow = T)  
 #Output matrix H (measurement equation): 
 h = matrix(c(1,0,0,1,0,0,0,
              0,0,0,0,0,0,1,
              0,0,0,0,1,1,0), nrow = 3, byrow = T)
 #Output measurement noise distribution matrix R:
 r = diag(0, 3,3)
 
 
 #Data for the model. Input: yt = (log gdp, inf, real rate); output = (rer_depreciation, import_inflation)
 data_dse_k = TSdata(input = data_kalman_ts[1:eval(nrow(data_kalman_ts)-1),c(3,2)], output = data_kalman_ts[2:nrow(data_kalman_ts),c(4,1,5)])
 windows()
 tfplot(data_dse_k)

 #Create grid for initial state values: 
 grid_gdp = seq(quantile(data_kalman$lgdp_tc, 0.25),quantile(data_kalman$lgdp_tc, 0.75), length.out = 3)
 grid_g = seq(0.01, 0.05, length.out = 3)
 grid_g_p = seq(0.01, 0.05, length.out = 3)
 grid_gap = seq(-0.05, 0.05, length.out = 3)
 grid_r_trend = seq(0.005, 0.04, length.out = 3)
 grid_r_gap = seq(-0.03, 0.03, length.out = 3)
 grid_inflation = seq(0.03, 0.1, length.out = 3)
 
 combinations = as.data.frame(crossing(gdp = 1:3, g = 1:3, g_p = 1:3, gap = 1:3, r_trend = 1:3, r_gap = 1:3, inflation = 1:3))
 dim(combinations)
 rmse_kalman = NULL
 negloglik_kalman = NULL
 ptm <- proc.time()
 for(i in 1:nrow(combinations)){
   #Construct state-space model (no innovations form):
      z0_vector = c(grid_gdp[combinations[1,i]],grid_g[combinations[2,i]],grid_g_p[combinations[3,i]],grid_gap[combinations[4,i]],grid_r_trend[combinations[5,i]]
                 ,grid_r_gap[combinations[6,i]],grid_inflation[combinations[7,i]])
   model_k1= SS(F.= f, G = g, H = h, Q = q, R = r, z0 = z0_vector) #c(11.5, 0.03, 0.03, 0, 0.015, 0, 0.03))
   
   #Estimation of model by maximum likelihood:
   model_k1_ml = estMaxLik(model_k1, data_dse_k, algorithm = "optim",
                           algorithm.args=list(method="L-BFGS-B", upper=c(0.99,0.99,0.99,0.99,0,0.99,Inf,0.99), lower=c(0,0,0,0,-Inf,0,0,0), hessian=TRUE,
                                               control=list(maxit=10000))) 
   summary_model_k1_ml = summary(model_k1_ml)
   loglik_aux = summary_model_k1_ml$estimates$l 
   rmse_aux = summary_model_k1_ml$estimates$rmse
   
   negloglik_kalman = rbind(negloglik_kalman, loglik_aux)
   rmse_kalman = rbind(rmse_kalman, rmse_aux)
 }
 proc.time() - ptm 
 combinations[7,]
 i = 7
 z0_vector = c(grid_gdp[combinations[1,i]],grid_g[combinations[2,i]],grid_g_p[combinations[3,i]],grid_gap[combinations[4,i]],grid_r_trend[combinations[5,i]]
               ,grid_r_gap[combinations[6,i]],grid_inflation[combinations[7,i]])
 #Selected z0 = ()
 ptm <- proc.time()
 #Construct state-space model (no innovations form):
 model_k1= SS(F.= f, G = g, H = h, Q = q, R = r, z0 = z0_vector) #c(11.5, 0.03, 0.03, 0, 0.015, 0, 0.03))
  
 #Estimation of model by maximum likelihood:
 model_k1_ml = estMaxLik(model_k1, data_dse_k, algorithm = "optim",
                         algorithm.args=list(method="L-BFGS-B", upper=c(0.99,0.99,0.99,0.99,0,0.99,Inf,0.99), lower=c(0,0,0,0,-Inf,0,0,0), hessian=TRUE,
                                             control=list(maxit=10000))) 
 
 summary(model_k1_ml)
 summary_model_k1_ml = summary(model_k1_ml)
 summary_model_k1_ml$estimates$l
 print(model_k1_ml)
 proc.time() - ptm 
 
 #Smoother for obtaining the state filtered values: 
 smoother_model1 = smoother(model_k1_ml)
 print(smoother_model1$estimates)
 print(smoother_model1$smooth$state)
 windows()
 plot(smoother_model1$estimates$pred[,1], type = "l")
 lines(data_dse_k$output[,1])
 plot(smoother_model1$smooth$state[,1], type = "l")
 #Potential output and output GAP:
 y_p = smoother_model1$smooth$state[,1]
 y_gap= smoother_model1$smooth$state[,4]
 windows()
 plot(y_gap, type = 'l')
 lines(data_dse_k$output[,1], type = 'l', col = 'blue')
 plot(smoother_model1$smooth$state[,1])
 
 write.csv(as.data.frame(cbind(y_p,y_gap)), "producto_potencial_kf.csv")
 windows()
 tfplot(state(model_k1_ml))

# ------------------- Kalman Filter with the DLM package --------------------- # 
 
 #Define function for MLE estimation: 
 my_dlmfc <- function(par = c(gamma_0, beta_1, beta_2, alpha_1, alpha_2, sigma_yp, sigma_g, sigma_gb, sigma_y, sigma_rb, sigma_zr, sigma_pi)){
   gamma_0 = par[1]
   beta_1 = par[2]
   beta_2 = par[3]
   alpha_1 = par[4]
   alpha_2 = par[5]
   sigma_yp = par[6]
   sigma_g = par[7]
   sigma_gb = par[8]
   sigma_y = par[9]
   sigma_rb = par[10]
   sigma_zr = par[11]
   sigma_pi = par[12]
   
   #Transition matrix G
   GG =  matrix(c(1,1,0,0,0,0,0,
                  0, gamma_0, 1-gamma_0,0,0,0,0,
                  0,0,1,0,0,0,0,
                  0,0,0,beta_1,0,-beta_2,0,
                  0,0,0,0,1,0,0,
                  0,0,0,0,0,0,0,
                  0,0,0,alpha_2,0,0,alpha_1), nrow = 7, byrow = T)
   
   #Measurement equation matrix
   FF =   matrix(c(1,0,0,1,0,0,0,
                   0,0,0,0,0,0,1,
                   0,0,0,0,1,1,0), nrow = 3, byrow = T)
   
   #Measurement equation var-cov matrix: 
   V = diag(0, 3,3)
   
   #State equation var-cov matrix: 
   W  = matrix(c(sigma_yp,0,0,0,0,0,0,
                  0,sigma_g,0,0,0,0,0,
                  0,0,sigma_gb,0,0,0,0,
                  0,0,0,sigma_y,0,0,0,
                  0,0,0,0,sigma_rb,0,0,
                  0,0,0,0,0,sigma_zr,0,
                  0,0,0,0,0,0,sigma_pi), nrow = 7, byrow = T)
   m0 = matrix(c(0,0,0,0,0,0,0), nrow = 7)
   C0 = matrix(rep(0,49), nrow = 7)
   my_dlm1 = dlm(FF=FF,V=V,GG=GG,W=exp(W),V=V,m0=m0,C0=C0)
 }

 #Data for estimation (GDP, core inflation, real interest rate): 
 data_dlm = cbind(data$gdp, data$inflation_cpi_core, data$dtf_real)
 
 my_dlm_mle = dlmMLE(data_dlm, parm = c(0.5, 0.9, 0.4,0.7,0.5,1,1,1,1,1,1,1),lower=c(0,0,0,0,0,-Inf, -Inf, -Inf, -Inf, -Inf, -Inf,-Inf), my_dlmfc)
 my_dlm_mle$convergence
 print(my_dlm_mle)
 model.fit = my_dlmfc(my_dlm_mle$par)
 print(model.fit)
 my_dlm_smoothed <- dlmSmooth(data_dlm, model.fit) 
 my_dlm_smoothed <- my_dlm_smoothed$s
 #Output Gap 
 output_gap_dlm = my_dlm_smoothed[,4]
 y_p_dlm = my_dlm_smoothed[,1]
 plot(y_p_dlm)
 
 
  

#
 
 
 
 
 
 
 
 
 
 





