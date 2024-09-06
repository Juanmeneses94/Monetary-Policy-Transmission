# ---------------------------------------------------------------------------- #
# Thesis: "Nonlinearity, State-Dependency, and Asymmetries in Monetary Policy
# Transmission Mechanisms"
# Msc in Economics, Universidad de los Andes
#
# Methodologies included in the code are: 
# - Non Linear Least Squares Estimate of Monetary Transmission Mechanism System 
# - Testing for Nonlinearity and State Dependency (LM Type-Tests)
# - Selection of transition function
# - Estimates of STR models and diagnosis tests
# - Estimates of GIRF (densities)
#
# - Author: Juan Camilo Meneses Cortés - Universidad de los Andes
# - Version: 1.0 
# - Last update: 22-01-2021
# ---------------------------------------------------------------------------- #

#Libraries 
library(readxl)
library(nlme)
library(nlstools)
library(regtools)
library(tseries)
library(lmtest)
library(astsa)
library(xts)
library(stats4)
library(aod)
require(devtools)      
library(numDeriv)
library(tidyr)
library(Hmisc)
library(nlmrt)
library(NonlinearTSA)
library(urca)
library(sandwich)
library(FinTS)
library(minpack.lm)
library(matlib)
library(parallel)
library(doParallel)
library(doSNOW)


#Working Directory
setwd("C:/Users/juanm/OneDrive - Universidad de los Andes/Juan Camilo Meneses/Personal/Tesis PEG/Datos")

#Read Data: 
 data = read_excel("Datos.xlsx", sheet = "Data", col_names = T, skip = 0)
 names(data)
  # ts data type 
  data_ts = ts(data[,-1], start = c(1994,1), frequency = 4)
  # xts data type
  data_xts = xts(data[,-1], as.Date(data$Fecha))

# ------------------------ Unit Root Tests ------------------------------------ 
 
 #Standard Unit Root Test  
    #Endogenous Model Variables: 
       #Augmented Dickey - Fuller Test: 
       summary(ur.df(na.omit(data_xts$output_gap), type = "drift", selectlags = "AIC"))
       summary(ur.df(na.omit(data_xts$inflation_core_15), type = "drift", selectlags = "AIC"))
       summary(ur.df(na.omit(data_xts$repo), type = "drift", selectlags = "AIC"))
       #ERS Test: 
       summary(ur.ers(na.omit(data_xts$output_gap), type = "P-test", model = "constant", lag.max = 8))
       summary(ur.ers(na.omit(data_xts$inflation_core_15), type = "P-test", model = "constant", lag.max = 8))
       summary(ur.ers(na.omit(data_xts$repo), type = "P-test", model = "constant", lag.max = 8))
       #KPSS Test:
       summary(ur.kpss(na.omit(data_xts$output_gap), type = "mu", lags = "short"))
       summary(ur.kpss(na.omit(data_xts$inflation_core_15), type = "mu", lags = "short"))
       summary(ur.kpss(na.omit(data_xts$repo), type = "mu", lags = "short"))
       
       
    #Exogenous Model Variables: 
       #ADF Test: 
       summary(ur.df(na.omit(data_xts$import_inflation_ppi), type = "drift", selectlags = "AIC"))
       summary(ur.df(na.omit(data_xts$rer_gap), type = "drift", selectlags = "AIC"))
       summary(ur.df(na.omit(data_xts$real_interest_gap), type = "drift", selectlags = "AIC"))
       summary(ur.df(na.omit(data_xts$brent_real_gap), type = "drift", selectlags = "AIC"))
       summary(ur.df(na.omit(data_xts$output_gap_us), type = "drift", selectlags = "AIC"))
       summary(ur.df(na.omit(data_xts$inflation_gap), type = "drift", selectlags = "AIC"))
       #ERS Test: 
       summary(ur.ers(na.omit(data_xts$import_inflation_ppi), type = "P-test", model = "constant", lag.max = 8))
       summary(ur.ers(na.omit(data_xts$rer_gap), type = "P-test", model = "constant", lag.max = 8))
       summary(ur.ers(na.omit(data_xts$real_interest_gap), type = "P-test", model = "constant", lag.max = 8))
       summary(ur.ers(na.omit(data_xts$brent_real_gap), type = "P-test", model = "constant", lag.max = 8))
       summary(ur.ers(na.omit(data_xts$output_gap_us), type = "P-test", model = "constant", lag.max = 8))
       #KPSS Test: 
       summary(ur.kpss(na.omit(data_xts$import_inflation_ppi), type = "mu", lags = "short"))
       summary(ur.kpss(na.omit(data_xts$rer_gap), type = "mu", lags = "short"))
       summary(ur.kpss(na.omit(data_xts$real_interest_gap), type = "mu", lags = "short"))
       summary(ur.kpss(na.omit(data_xts$brent_real_gap), type = "mu", lags = "short"))
       summary(ur.kpss(na.omit(data_xts$output_gap_us), type = "mu", lags = "short"))
    
 #Enders and Ludlow (1999-2000) Non-linear decay test with Fourier
    enders_ludlow = function(series, p){
       demeaned_series = residuals(lm(na.omit(series)~1))
       diff_demeaned_series = diff(demeaned_series)
       base_el = cbind(diff_demeaned_series, lag(demeaned_series))
       final_el_test = NULL
       for(w in 1:p){
          for(i in 1:w){base_el = cbind(base_el,lag(base_el[,1],i))}
          base_el = na.omit(base_el)
          T = nrow(base_el)
          time_trend = seq(1,T,1)
          k_grid = seq(1,T/2,1)
          ssr_el = NULL
          for(k in 1:length(k_grid)){
          #Construction of the series of alpha():   
             sin_variable = sin(((2*pi*k)/T)*time_trend)*base_el[,2]
             cos_variable = cos(((2*pi*k)/T)*time_trend)*base_el[,2]
          #OLS regression without intercept
          ols_est = lm(base_el[,1]~ base_el[,2]+sin_variable + cos_variable + base_el[,3:eval(3+w-1)]-1, data = base_el)    
          ssr_el_aux <- sum((predict(ols_est) - base_el[,1])^2)  ## residual sum of squares
          ssr_el = rbind(ssr_el, ssr_el_aux)
          }
          #Selection of k that minimizes ssr
          k_selected = k_grid[which(ssr_el == min(ssr_el))]
          #Estimation of final model
          sin_variable_final = sin(((2*pi*k_selected)/T)*time_trend)*base_el[,2]
          cos_variable_final = cos(((2*pi*k_selected)/T)*time_trend)*base_el[,2]
          #OLS regression without intercept
          ols_est_final = lm(base_el[,1]~ base_el[,2]+sin_variable_final + cos_variable_final + base_el[,3:eval(3+w-1)]-1, data = base_el)    
          ssr_el_final <- sum((predict(ols_est_final) - base_el[,1])^2)  ## residual sum of squares
          summary_ols = summary(ols_est_final)
          var_cov_matrix_final = as.data.frame(vcov(ols_est_final))
          aic_model = AIC(ols_est_final)
          bic_model = BIC(ols_est_final)
          lb_test = Box.test(residuals(ols_est_final), lag = T/4)$p.value
            #C test: 
            C_test = summary_ols$coefficients[1,3]
            
            #CR test 
            c = ols_est_final$coefficients[1]
            a1 = ols_est_final$coefficients[2]
            a1_se = summary_ols$coefficients[2,2]
            b1 = ols_est_final$coefficients[3]
            b1_se = summary_ols$coefficients[3,2]
            r = sqrt(a1^2 + b1^2)
            CR = (r^2)/4
            
            dr_da1 = a1*((a1^2 + b1^2)^(-0.5))
            dr_db1 = b1*((a1^2 + b1^2)^(-0.5))
            cov_a1_b1 = var_cov_matrix_final[3,3]
            est_asy_var_cr = (dr_da1^2)*(a1_se^2) + (dr_db1^2)*(b1_se^2) + 2*dr_da1*dr_db1*cov_a1_b1 
            se_cr = sqrt(est_asy_var_cr)
            
            CR_test = (CR - c)/se_cr 
            
            #F_all test: 
            ols_est_F_all = lm(base_el[,1]~ base_el[,3:eval(3+w-1)]-1, data = base_el)
            ssr_el_F_all <- sum((predict(ols_est_F_all) - base_el[,1])^2)  ## residual sum of squares
            F_test_F_all = ((ssr_el_F_all-ssr_el_final)/3)/(ssr_el_final/(nrow(base_el)-length(ols_est_final$coefficients)))   
   
            #F_trig test: 
            ols_est_F_trig = lm(base_el[,1]~ base_el[,2]+ base_el[,3:eval(3+w-1)]-1, data = base_el)
            ssr_el_F_trig <- sum((predict(ols_est_F_trig) - base_el[,1])^2)  ## residual sum of squares
            F_test_F_trig = ((ssr_el_F_trig-ssr_el_final)/2)/(ssr_el_final/(nrow(base_el)-length(ols_est_final$coefficients)))   
   
            final_el_aux = cbind(F_test_F_all, F_test_F_trig, C_test, CR_test, aic_model, bic_model, lb_test)
            final_el_test = rbind(final_el_test, final_el_aux)
            }
            p_aic = which(final_el_test[,5]==min(final_el_test[,5]))
            p_bic = which(final_el_test[,6]==min(final_el_test[,6]))
            final_el_test_aic = rbind(as.data.frame(final_el_test[p_aic,]), p_aic)
            final_el_test_bic = rbind(as.data.frame(final_el_test[p_bic,]), p_bic)
            
            salida_final = cbind(final_el_test_aic,final_el_test_bic)
            rownames(salida_final) = c("F_all", "F_trig", "C", "CR", "AIC", "BIC", "LB", "p")
            colnames(salida_final) = c("AIC", "BIC")
            return(salida_final)
    }
       #Niveles
       enders_ludlow(data_xts$output_gap["1999-03-01/"], 8)
       enders_ludlow(data_xts$inflation_core_15["2002-03-01/"], 8)
       enders_ludlow(data_xts$repo["2002-03-01/"],8)
       #Primera Diferencia
       enders_ludlow(diff(data_xts$output_gap["1999-03-01/"]), 8)
       enders_ludlow(diff(data_xts$inflation_core_15["1999-03-01/"]),8)       
       enders_ludlow(diff(data_xts$repo["1999-03-01/"]),8)           
 
      
 # Tests based on package Nonlinear TSA
 Sollis_2004_unit_root(na.omit(data_xts$inflation_core_15),1, max_lags = 8)
 Cuestas_Garratt_unit_root(na.omit(data_xts$inflation_core_15), max_lags = 8, 1)
 Hu_Chen_Unit_Root(na.omit(data_xts$inflation_core_15),1,lags = 8, 1)
   
# ------------------------- Least Square Estimates -----------------------------
data_lm = cbind.xts(data_xts$output_gap, data_xts$inflation_core_15,
                data_xts$repo, lag(data_xts$output_gap,1), 
                lag(data_xts$inflation_core_15,1), lag(data_xts$inflation_core_15,2),lag(data_xts$repo,1),lag(data_xts$repo,6),data_xts$rer_gap,
                lag(data_xts$rer_gap,1), lag(data_xts$import_inflation_ppi,2),
                data_xts$real_interest_gap,lag(data_xts$real_interest_gap),lag(data_xts$real_interest_gap,2),
                lag(data_xts$real_interest_gap,3),lag(data_xts$real_interest_gap,4),lag(data_xts$real_interest_gap,5),
                lag(data_xts$real_interest_gap,6), data_xts$output_gap_us, data_xts$inflation_gap, data_xts$brent_real_gap)  
colnames(data_lm) = c("output_gap", "inflation", "interest_rate","output_gap.l1","inflation.l1", "inflation.l2",
                       "interest_rate.l1","interest_rate.l6","rer_gap", "rer_gap.l1", "inflation_ext.l2",
                      "real_interest_gap","real_interest_gap.l1","real_interest_gap.l2","real_interest_gap.l3","real_interest_gap.l4",
                      "real_interest_gap.l5", "real_interest_gap.l6", "us_output_gap", "inflation_gap", "brent_gap")
data_lm = na.omit(data_lm["2002-09-01/2019-12-01"])
nrow(data_lm)

head(data_lm)
tail(data_lm)
plot(lag(data_xts$real_interest_gap,4))
lines(data_xts$output_gap)
plot(data_xts$repo)
cor(data_xts$output_gap["2002-09-01/2019-12-01"],lag(data_xts$real_interest_gap,4)["2002-09-01/2019-12-01"], method = "pearson")

#Inflation Equation
  inflation_eq_lm = lm(inflation ~ 1+inflation.l1+output_gap+output_gap.l1+real_interest_gap.l2+
                     +rer_gap+inflation_ext.l2, data = data_lm)
   summary(inflation_eq_lm)
   AIC_inf_lm = AIC(inflation_eq_lm)
   BIC_inf_lm = BIC(inflation_eq_lm)
    #Test for the model: 
    #Normality test
    jarque.bera.test(residuals(inflation_eq_lm)) #Normal errors
    # Heterocedasticity test. 
    bptest(inflation_eq_lm) #Homoskedasticity at 5%
    #Serial Correlation of residuals
    bgtest(inflation_eq_lm, order = 12) #Serial correlation
    acf2(residuals(inflation_eq_lm))
    adjreg = sarima(residuals(inflation_eq_lm), 1,0,0, xreg = NULL) # This is the adjustment regression with MA(1) residuals
    adjreg # Results of adjustment regression. White noise should be suggested.
    acf2(resid(adjreg$fit), main = "Análisis de residuales estimación OLS ecuación inflación")
    
    
 #Output gap Equation: 
    outputgap_eq_lm = lm(output_gap ~ 1+output_gap.l1
                         +real_interest_gap.l2+rer_gap+us_output_gap + brent_gap, data = data_lm)
    summary(outputgap_eq_lm)
    AIC_outputgap_lm = AIC(outputgap_eq_lm)
    BIC_outputgap_lm = BIC(outputgap_eq_lm)
    #Test for the model: 
    #Normality test
    jarque.bera.test(residuals(outputgap_eq_lm)) #Normal errors at 10%
    # Heterocedasticity test. 
    bptest(outputgap_eq_lm) #Suggest Homoskedasticity
    #Serial Correlation of residuals
    acf2(residuals(outputgap_eq_lm))
    bgtest(outputgap_eq_lm, order = 12) #Suggest Serial Correlation
    adjreg = sarima(residuals(outputgap_eq_lm), 1,0,0, xreg = NULL) # This is the adjustment regression with AR(1) residuals
    adjreg # Results of adjustment regression. White noise should be suggested.
    acf2(resid(adjreg$fit), main = "Residuales OLS ecuación Brecha PIB")
    
  #Interest Rate Equation (Policy Rule): 
    interest_eq_lm = lm(interest_rate ~ 1+inflation_gap+output_gap+interest_rate.l1, data = data_lm)
    summary(interest_eq_lm)
    AIC_interest_lm = AIC(interest_eq_lm)
    BIC_interest_lm = BIC(interest_eq_lm)
    #Test for the model: 
    #Normality test
    jarque.bera.test(residuals(interest_eq_lm)) #No normal errors
    # Heterocedasticity test. 
    bptest(interest_eq_lm) #Suggest Heteroskedasticity 
    #Serial Correlation of residuals
    bgtest(interest_eq_lm, order = 12) # Suggest serial correlation
    acf2(residuals(interest_eq_lm))
    adjreg = sarima(residuals(interest_eq_lm), 1,0,0, xreg = NULL) # This is the adjustment regression with AR(1) residuals
    adjreg # Results of adjustment regression. White noise should be suggested.
    acf2(resid(adjreg$fit), main = "Residuales OLS ecuación Tasa Repo")
  
# -------------------- Nonlinear Least Square Estimates ------------------------
data_nls = cbind.xts(data_xts$output_gap, data_xts$inflation_core_15,
                     data_xts$repo, lag(data_xts$output_gap,1), 
                     lag(data_xts$inflation_core_15,1), lag(data_xts$inflation_core_15,2),lag(data_xts$repo,1),lag(data_xts$repo,2),data_xts$rer_gap,
                     lag(data_xts$rer_gap,1), lag(data_xts$import_inflation_ppi,2),
                     data_xts$real_interest_gap,lag(data_xts$real_interest_gap),lag(data_xts$real_interest_gap,2),
                     lag(data_xts$real_interest_gap,3),lag(data_xts$real_interest_gap,4),lag(data_xts$real_interest_gap,5),
                     lag(data_xts$real_interest_gap,6), data_xts$output_gap_us, data_xts$inflation_gap, data_xts$brent_real_gap)    
colnames(data_nls) = c("output_gap", "inflation", "interest_rate","output_gap.l1","inflation.l1", "inflation.l2",
                       "interest_rate.l1","interest_rate.l2","rer_gap", "rer_gap.l1", "inflation_ext.l2",
                       "real_interest_gap","real_interest_gap.l1","real_interest_gap.l2","real_interest_gap.l3","real_interest_gap.l4",
                       "real_interest_gap.l5", "real_interest_gap.l6", "us_output_gap", "inflation_gap", "brent_gap")
data_nls = na.omit(data_nls["2002-09-01/2019-12-01"])
obs = nrow(data_nls)
head(data_nls)
tail(data_nls)
    
data_nls = as.data.frame(data_nls)

#Inflation Equation
acf2(data_nls$inflation)

inflation_eq_nls=try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4
                        +b5*x5+b6*x6
                        -rho*b1*x7-rho*b2*x8-rho*b3*x9-rho*b4*x10
                        -rho*b5*x11-rho*b6*x12, 
                        start = list(rho = 0.01, b0 = 0.01, b1 = 0.67, b2 = -0.11, b3= 0.38, b4= -0.09, b5= 0.084, b6 = 0.02),
                        data = list(y = data_nls$inflation[2:obs], y_l = data_nls$inflation[1:obs-1], 
                                     x1 = data_nls$inflation.l1[2:obs], x2 = data_nls$output_gap[2:obs], x3 = data_nls$output_gap.l1[2:obs], 
                                     x4 = data_nls$real_interest_gap.l2[2:obs], x5 = data_nls$rer_gap[2:obs], x6 = data_nls$inflation_ext.l2[2:obs],
                                     x7 = data_nls$inflation.l1[1:obs-1], x8 = data_nls$output_gap[1:obs-1], x9 = data_nls$output_gap.l1[1:obs-1], 
                                     x10 = data_nls$real_interest_gap.l2[1:obs-1], x11 = data_nls$rer_gap[1:obs-1], x12 = data_nls$inflation_ext.l2[1:obs-1]
                                     ), 
                        trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200)))
print(summary(inflation_eq_nls))   
overview(inflation_eq_nls)
plot(inflation_eq_nls)
residuals_inf_nls = nlsResiduals(inflation_eq_nls)
plot(residuals_inf_nls)
test.nlsResiduals(residuals_inf_nls)
predicted_inf_nls = predict(inflation_eq_nls)
plot(predicted_inf_nls, type = 'l')
lines(data_nls$inflation[2:obs], col = 'blue')

   #Test for the model 
      #Normality test
      jarque.bera.test(residuals(inflation_eq_nls)) #Normality   
      shapiro.test(residuals(inflation_eq_nls))
      
      #Independence or residuals 
      LB_inf_nls = NULL
      for(i in 1:16){
           LB = Box.test(residuals(inflation_eq_nls), lag = i, type = "Ljung-Box")$p.value #Independence
           LB_inf_nls = rbind(LB_inf_nls,LB)
      }
      print(LB_inf_nls)
      # Heterocedasticity test (manual): 
      nlshc(inflation_eq_nls)
      summary(lm(residuals_inf_nls$resi1[,2]^2 ~ 1 +data_nls$inflation.l1[2:obs]+data_nls$real_interest_gap.l2[2:obs]+
                    data_nls$output_gap[2:obs]+data_nls$output_gap.l1[2:obs]
                 +data_nls$rer_gap[2:obs]+data_nls$inflation_ext.l2[2:obs]))
      # Serial Correlation: 
      lm_aux_bp = lm(residuals_inf_nls$resi1[,2]~ 1)
      bgtest(lm_aux_bp, order = 12)
      acf2(residuals_inf_nls$resi1[,2])
   
   ssr_1_inf_nls <- sum((predicted_inf_nls - data_nls$inflation[2:obs]) ^ 2)  ## residual sum of squares
   regss_inf_nls <- sum((predicted_inf_nls - mean(predicted_inf_nls)) ^ 2) ## regression sum of squares
   tss_inf_nls <- sum((data_nls$inflation[2:obs] - mean(data_nls$inflation[2:obs])) ^ 2)  ## total sum of squares
   r_squared_inf_nls = regss_inf_nls / tss_inf_nls
   print(r_squared_inf_nls)
   BIC_inf_nls = BIC(inflation_eq_nls)
   AIC_inf_nls = AIC(inflation_eq_nls)
   print(BIC_inf_nls)
   print(AIC_inf_nls)
  
#Output Gap Equation
acf2(data_nls$output_gap)

   outputgap_eq_nls=try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4+b5*x5
                            -rho*b1*x6-rho*b2*x7-rho*b3*x8-rho*b4*x9-rho*b5*x10,
                            start = list(rho = 0.03, b0 = -0.0015, b1 = 0.88, b2 = -0.023, b3= -0.084, b4= 0.23, b5 = 0.004),
                            data = list(y = data_nls$output_gap[2:obs], y_l = data_nls$output_gap[1:obs-1], x1 = data_nls$output_gap.l1[2:obs],
                                        x2 = data_nls$real_interest_gap.l2[2:obs], x3 = data_nls$us_output_gap[2:obs], x4 = data_nls$rer_gap[2:obs], x5 = data_nls$brent_gap[2:obs], 
                                        x6 = data_nls$output_gap.l1[1:obs-1], x7 = data_nls$real_interest_gap.l2[1:obs-1], x8 = data_nls$us_output_gap[1:obs-1], x9 = data_nls$rer_gap[1:obs-1], x10 = data_nls$brent_gap[1:obs-1]
                                        )
                            ,trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200)))
   print(summary(outputgap_eq_nls))   
   overview(outputgap_eq_nls)
   plot(outputgap_eq_nls)
   residuals_outputgap_nls = nlsResiduals(outputgap_eq_nls)
   plot(residuals_outputgap_nls)
   test.nlsResiduals(residuals_outputgap_nls)
   predicted_outputgap_nls = predict(outputgap_eq_nls)
   plot(predicted_outputgap_nls, type = 'l', col = "blue")
   lines(data_nls$output_gap[2:obs])
   residuals(outputgap_eq_nls)
   #Test for the model 
      #Normality test
      jarque.bera.test(residuals(outputgap_eq_nls)) #Normality   
      shapiro.test(residuals(outputgap_eq_nls))
      #Independence or residuals 
      LB_output_nls = NULL
      for(i in 1:16){
         LB = Box.test(residuals(outputgap_eq_nls), lag = i, type = "Ljung-Box")$p.value #Independence
         LB_output_nls = rbind(LB_output_nls,LB)
      }
      print(LB_output_nls)
      # Heterocedasticity test (manual): 
      nlshc(outputgap_eq_nls)
      summary(lm(residuals_outputgap_nls$resi1[,2]^2 ~ 1+data_nls$output_gap.l1[2:obs]+data_nls$real_interest_gap.l2[2:obs]+
                    data_nls$rer_gap[2:obs]+data_nls$us_output_gap[2:obs]+data_nls$brent_gap[2:obs])) # Homoskedasticiy at 5%
      # Serial Correlation: 
      lm_aux_bp = lm(residuals_outputgap_nls$resi1[,2]~ 1)
      bgtest(lm_aux_bp, order = 12)
      acf2(residuals_outputgap_nls$resi1[,2])
      
   ssr_1_outputgap_nls <- sum((predicted_outputgap_nls - data_nls$output_gap[2:obs]) ^ 2)  ## residual sum of squares
   regss_outputgap_nls <- sum((predicted_outputgap_nls - mean(predicted_outputgap_nls)) ^ 2) ## regression sum of squares
   tss_outputgap_nls <- sum((data_nls$output_gap[2:obs] - mean(data_nls$output_gap[2:obs])) ^ 2)  ## total sum of squares
   r_squared_outputgap_nls = regss_outputgap_nls/tss_outputgap_nls
   print(r_squared_outputgap_nls)
   BIC_outputgap_nls = BIC(outputgap_eq_nls)
   AIC_outputgap_nls = AIC(outputgap_eq_nls)
   print(BIC_outputgap_nls)
   print(AIC_outputgap_nls)   
   
#Interest Rate Equation: 
acf2(data_nls$interest_rate)
plot(data_nls$interest_rate, type = "l")
   
   interest_eq_nls=try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+
                               -rho*b1*x4-rho*b2*x5-rho*b3*x6,
                            start = list(rho = 0.03, b0 = 0.0, b1 = 0.9, b2 = 0.4, b3= 0.15),
                            data = list(y = data_nls$interest_rate[2:obs], y_l = data_nls$interest_rate[1:obs-1], x1 = data_nls$interest_rate.l1[2:obs],
                                        x2 = data_nls$inflation_gap[2:obs], x3 = data_nls$output_gap[2:obs], x4 = data_nls$interest_rate.l1[1:obs-1], 
                                        x5 = data_nls$inflation_gap[1:obs-1], x6 = data_nls$output_gap[1:obs-1])
                            ,trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200)))
   print(summary(interest_eq_nls))   
   overview(interest_eq_nls)
   plot(interest_eq_nls)
   residuals_interest_nls = nlsResiduals(interest_eq_nls)
   plot(residuals_interest_nls)
   test.nlsResiduals(residuals_interest_nls)
   predicted_interest_nls = predict(interest_eq_nls)
   plot(predicted_interest_nls, type = 'l', col = "gray")
   lines(data_nls$interest_rate[2:obs])
   residuals(interest_eq_nls)
   #Test for the model 
   #Normality test
   jarque.bera.test(residuals(interest_eq_nls)) #Normality   
   shapiro.test(residuals(interest_eq_nls))
   #Independence or residuals 
   LB_interest_nls = NULL
   for(i in 1:16){
      LB = Box.test(residuals(interest_eq_nls), lag = i, type = "Ljung-Box")$p.value #Independence
      LB_interest_nls = rbind(LB_interest_nls,LB)
   }
   print(LB_interest_nls)
   # Heterocedasticity test (manual): 
   nlshc(interest_eq_nls)
   summary(lm(residuals_interest_nls$resi1[,2]^2 ~ 1 +data_nls$interest_rate.l1[2:obs]+data_nls$inflation_gap[2:obs]+
                 data_nls$output_gap[2:obs])) # Heteroskedasticity at 5%
   sandwich(interest_eq_nls)
   vcov(interest_eq_nls)
   coeftest(interest_eq_nls, vcov = sandwich(interest_eq_nls))
   # Serial Correlation: 
   lm_aux_bp = lm(residuals_interest_nls$resi1[,2]~ 1)
   bgtest(lm_aux_bp, order = 12)
   acf2(residuals_interest_nls$resi1[,2])
   
   ssr_1_interest_nls <- sum((predicted_interest_nls - data_nls$interest_rate[2:obs]) ^ 2)  ## residual sum of squares
   regss_interest_nls <- sum((predicted_interest_nls - mean(predicted_interest_nls)) ^ 2) ## regression sum of squares
   tss_interest_nls <- sum((data_nls$interest_rate[2:obs] - mean(data_nls$interest_rate[2:obs])) ^ 2)  ## total sum of squares
   r_squared_interest_nls = regss_interest_nls / tss_interest_nls
   print(r_squared_interest_nls)
   BIC_interest_nls = BIC(interest_eq_nls)
   AIC_interest_nls = AIC(interest_eq_nls)
   print(BIC_interest_nls)
   print(AIC_interest_nls)  
   
# ----------------------- Correlation Test of Residuals ------------------------

#Non- linear least square estimator
rcorr(as.matrix(cbind(residuals_inf_nls$resi1[,2], residuals_outputgap_nls$resi1[,2], residuals_interest_nls$resi1[,2])), type = "pearson")

# ------------------------ Direct Non-linearity Tests --------------------------        
# ------------------- Ramsey Test ------------------ 
# The test are performed over the NLS or GLS estimations
# Inflation Equation: 
data_inf_ramsey_nls = cbind(data_nls[2:obs,], predicted_inf_nls, predicted_inf_nls^2, predicted_inf_nls^3, predicted_inf_nls^4)
obs_ramsey = nrow(data_inf_ramsey_nls)
names(data_inf_ramsey_nls)
colnames(data_inf_ramsey_nls) =  c("output_gap", "inflation", "interest_rate","output_gap.l1","inflation.l1", "inflation.l2",
                               "interest_rate.l1","interest_rate.l2","rer_gap", "rer_gap.l1", "inflation_ext.l2",
                               "real_interest_gap","real_interest_gap.l1","real_interest_gap.l2","real_interest_gap.l3" ,"real_interest_gap.l4","real_interest_gap.l5", "real_interest_gap.l6",
                               "us_output_gap","inflation_gap","brent_gap",
                               "predicted_inf", "predicted_inf.2", "predicted_inf.3", "predicted_inf.4") 
  
  
  reset_inf = resettest(inflation ~ 1+inflation.l1+output_gap+output_gap.l1+real_interest_gap.l2
                                   +rer_gap+inflation_ext.l2
                        ,data = data_inf_ramsey_nls, type = "fitted", power = 1:4)
  reset_inf$p.value
  reset_regression_inf_nls = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4
                                     +b5*x5+b6*x6
                                     +b7*x7+b8*x8+b9*x9+b10*x10
                                     -rho*b1*x11-rho*b2*x12-rho*b3*x13-rho*b4*x14
                                     -rho*b5*x15-rho*b6*x16
                                     -rho*b7*x17-rho*b8*x18-rho*b9*x19-rho*b10*x10, 
                                     start = list(rho = 0.01, b0 = 0.01, b1 = 0.67, b2 = -0.11, b3= 0.38, b4= -0.09, b5= 0.084, b6 = 0.02, b7 = 0, b8 = 0, b9 = 0, b10 = 0),
                                     data = list(y = data_inf_ramsey_nls$inflation[2:obs_ramsey], y_l = data_inf_ramsey_nls$inflation[1:obs_ramsey-1], 
                                                 x1 = data_inf_ramsey_nls$inflation.l1[2:obs_ramsey], x2 = data_inf_ramsey_nls$output_gap[2:obs_ramsey], x3 = data_nls$output_gap.l1[2:obs_ramsey], 
                                                 x4 = data_inf_ramsey_nls$real_interest_gap.l2[2:obs_ramsey], x5 = data_inf_ramsey_nls$rer_gap[2:obs_ramsey], x6 = data_inf_ramsey_nls$inflation_ext.l2[2:obs_ramsey],
                                                 x7 = data_inf_ramsey_nls$predicted_inf[2:obs_ramsey],x8 = data_inf_ramsey_nls$predicted_inf.2[2:obs_ramsey], x9 = data_inf_ramsey_nls$predicted_inf.3[2:obs_ramsey],x10 = data_inf_ramsey_nls$predicted_inf.4[2:obs_ramsey],
                                                 x11 = data_nls$inflation.l1[1:obs_ramsey-1], x12 = data_nls$output_gap[1:obs_ramsey-1], x13 = data_nls$output_gap.l1[1:obs_ramsey-1], 
                                                 x14 = data_nls$real_interest_gap.l2[1:obs_ramsey-1], x15 = data_nls$rer_gap[1:obs_ramsey-1], x16 = data_nls$inflation_ext.l2[1:obs_ramsey-1],
                                                 x17 = data_inf_ramsey_nls$predicted_inf[1:obs_ramsey-1],x18 = data_inf_ramsey_nls$predicted_inf.2[1:obs_ramsey-1],x19 = data_inf_ramsey_nls$predicted_inf.3[1:obs_ramsey-1], x20 = data_inf_ramsey_nls$predicted_inf.4[1:obs_ramsey-1]
                                     ), 
                                     trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 100)))
  
  predicted_inf_nls_ramsey = predict(reset_regression_inf_nls)   
  ssr_inf_nls_ramsey <- sum((predicted_inf_nls_ramsey - data_inf_ramsey_nls$inflation[2:obs_ramsey]) ^ 2)  ## residual sum of squares
  #F-teset over with H0: predicted inf powers = 0
  F_inf_ramsey <-((ssr_1_inf_nls - ssr_inf_nls_ramsey)/(1+length(reset_regression_inf_nls$m$getPars())-length(inflation_eq_nls$m$getPars())))/(ssr_inf_nls_ramsey/(nrow(data_inf_ramsey_nls[2:obs_ramsey,])-length(reset_regression_inf_nls$m$getPars())))
  pval_F_inf_ramsey = pf(F_inf_ramsey,1+length(reset_regression_inf_nls$m$getPars())-length(inflation_eq_nls$m$getPars()),nrow(data_inf_ramsey_nls[2:obs_ramsey,])-length(reset_regression_inf_nls$m$getPars()), lower.tail = FALSE)
  #Chi-Square version of the RESET test: 
  Chi_inf_ramsey <-nrow(data_inf_ramsey_nls[2:obs_ramsey,])*((ssr_1_inf_nls - ssr_inf_nls_ramsey)/(ssr_1_inf_nls))
  print(Chi_inf_ramsey)
  pval_chi_inf_ramsey = pchisq(Chi_inf_ramsey,1+length(reset_regression_inf_nls$m$getPars())-length(inflation_eq_nls$m$getPars()), lower.tail = FALSE)
  
  #Test over RAMSEY Equation: 
  residuals_inf_ramsey = nlsResiduals(reset_regression_inf_nls)$resi1[,2]
    
    #Normality test
    jarque.bera.test(residuals(reset_regression_inf_nls)) #Normality   
    shapiro.test(residuals(reset_regression_inf_nls)) 
    # Heterocedasticity test (manual): 
    summary(lm(resid(reset_regression_inf_nls)^2 ~1+inflation.l1+output_gap+output_gap.l1+real_interest_gap.l2
               +rer_gap+inflation_ext.l2+predicted_inf+
                 predicted_inf.2+predicted_inf.3
               +predicted_inf.4, data = data_inf_ramsey_nls[2:obs_ramsey,])) #Homoskedasticiy at 1%
    # Serial Correlation: 
    lm_aux_bp = lm(resid(reset_regression_inf_nls)~ 1)
    bgtest(lm_aux_bp, order = 12)
    acf2(resid(reset_regression_inf_nls))
    
  
#Output Gap Equation:
  data_outputgap_ramsey_nls = cbind(data_nls[2:obs,], predicted_outputgap_nls, predicted_outputgap_nls^2, predicted_outputgap_nls^3, predicted_outputgap_nls^4)
  obs_ramsey = nrow(data_outputgap_ramsey_nls)
  names(data_outputgap_ramsey_nls)
  colnames(data_outputgap_ramsey_nls) =  c("output_gap", "inflation", "interest_rate","output_gap.l1","inflation.l1", "inflation.l2",
                                           "interest_rate.l1","interest_rate.l2","rer_gap", "rer_gap.l1", "inflation_ext.l2",
                                           "real_interest_gap","real_interest_gap.l1","real_interest_gap.l2","real_interest_gap.l3" ,"real_interest_gap.l4","real_interest_gap.l5", "real_interest_gap.l6",
                                           "us_output_gap","inflation_gap","brent_gap",
                                       "predicted_outputgap", "predicted_outputgap.2", "predicted_outputgap.3", "predicted_outputgap.4")  
  
  reset_regression_output_gap = resettest(output_gap ~ 1+output_gap.l1+real_interest_gap.l2+
                                            us_output_gap+rer_gap + brent_gap
                                          ,data = data_outputgap_ramsey_nls, type = "fitted", power = 1:4) 
  reset_regression_output_gap$p.value
  reset_regression_output_gap_nls = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4+b5*x5
                                            +b6*x6+b7*x7+b8*x8+b9*x9
                                            -rho*b1*x10-rho*b2*x11-rho*b3*x12-rho*b4*x13-rho*b5*x14
                                            -rho*b6*x15-rho*b7*x16-rho*b8*x17-rho*b9*x18,
                                            start = list(rho = 0.03, b0 = -0.0015, b1 = 0.88, b2 = -0.023, b3= -0.084, b4= 0.23, b5 = 0.004, b6 = 0, b7 = 0, b8 = 0, b9 = 0),
                                            data = list(y = data_outputgap_ramsey_nls$output_gap[2:obs_ramsey], y_l = data_outputgap_ramsey_nls$output_gap[1:obs_ramsey-1], x1 = data_outputgap_ramsey_nls$output_gap.l1[2:obs_ramsey],
                                                        x2 = data_outputgap_ramsey_nls$real_interest_gap.l2[2:obs_ramsey], x3 = data_outputgap_ramsey_nls$us_output_gap[2:obs_ramsey], x4 = data_outputgap_ramsey_nls$rer_gap[2:obs_ramsey], x5 = data_outputgap_ramsey_nls$brent_gap[2:obs_ramsey], 
                                                        x6 = data_outputgap_ramsey_nls$predicted_outputgap[2:obs_ramsey], x7 = data_outputgap_ramsey_nls$predicted_outputgap.2[2:obs_ramsey], x8 = data_outputgap_ramsey_nls$predicted_outputgap.3[2:obs_ramsey], x9 = data_outputgap_ramsey_nls$predicted_outputgap.4[2:obs_ramsey], 
                                                        x10 = data_outputgap_ramsey_nls$output_gap.l1[1:obs_ramsey-1], x11 = data_outputgap_ramsey_nls$real_interest_gap.l2[1:obs_ramsey-1], x12 = data_outputgap_ramsey_nls$us_output_gap[1:obs_ramsey-1], x13 = data_outputgap_ramsey_nls$rer_gap[1:obs_ramsey-1], x14 = data_outputgap_ramsey_nls$brent_gap[1:obs_ramsey-1], 
                                                        x15 = data_outputgap_ramsey_nls$predicted_outputgap[1:obs_ramsey-1], x16 = data_outputgap_ramsey_nls$predicted_outputgap.2[1:obs_ramsey-1], x17 = data_outputgap_ramsey_nls$predicted_outputgap.3[1:obs_ramsey-1], x18 = data_outputgap_ramsey_nls$predicted_outputgap.4[1:obs_ramsey-1]
                                            )
                                            ,trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200))) 
  
  summary(reset_regression_output_gap_nls)
  predicted_outputgap_nls_ramsey = predict(reset_regression_output_gap_nls)   
  ssr_outputgap_nls_ramsey <- sum((predicted_outputgap_nls_ramsey - data_outputgap_ramsey_nls$output_gap[2:obs_ramsey]) ^ 2)  ## residual sum of squares
  #F-teset over with H0: predicted inf powers = 0
  F_outputgap_ramsey <-((ssr_1_outputgap_nls - ssr_outputgap_nls_ramsey)/(1+length(reset_regression_output_gap_nls$m$getPars())-length(outputgap_eq_nls$m$getPars())))/(ssr_outputgap_nls_ramsey/(nrow(data_outputgap_ramsey_nls[2:obs_ramsey,])-length(reset_regression_output_gap_nls$m$getPars())))
  pf(F_outputgap_ramsey,1+length(reset_regression_output_gap_nls$m$getPars())-length(outputgap_eq_nls$m$getPars()),nrow(data_outputgap_ramsey_nls[2:obs_ramsey,])-length(reset_regression_output_gap_nls$m$getPars()), lower.tail = FALSE)
  #Chi-Square version of the RESET test: 
  Chi_outputgap_ramsey <-nrow(data_outputgap_ramsey_nls[2:obs_ramsey,])*((ssr_1_outputgap_nls - ssr_outputgap_nls_ramsey)/(ssr_1_outputgap_nls))
  print(Chi_outputgap_ramsey)
  pchisq(Chi_outputgap_ramsey,1+length(reset_regression_output_gap_nls$m$getPars())-length(outputgap_eq_nls$m$getPars()), lower.tail = FALSE)
  
  #Test over RAMSEY Equation: 
  residuals_outputgap_ramsey = nlsResiduals(reset_regression_output_gap_nls)$resi1[,2]
  
  #Normality test
  jarque.bera.test(residuals(reset_regression_output_gap_nls)) #Normality   
  shapiro.test(residuals(reset_regression_output_gap_nls)) 
  # Heterocedasticity test (manual): 
  summary(lm(resid(reset_regression_output_gap_nls)^2 ~1+output_gap.l1+real_interest_gap.l2+
             us_output_gap+rer_gap+brent_gap+
             predicted_outputgap+predicted_outputgap.2+predicted_outputgap.3
          +predicted_outputgap.4, data = data_outputgap_ramsey_nls[2:obs_ramsey,])) #Homoskedasticiy at 5%
  # Serial Correlation: 
  lm_aux_bp = lm(resid(reset_regression_output_gap_nls)~ 1)
  bgtest(lm_aux_bp, order = 12)
  acf2(resid(reset_regression_output_gap_nls)) 
  
    
#Interest Rate Equation:
  data_interest_ramsey_nls = cbind(data_nls[2:obs,], predicted_interest_nls, predicted_interest_nls^2, predicted_interest_nls^3, predicted_interest_nls^4)
  obs_ramsey = nrow(data_interest_ramsey_nls)
  names(data_interest_ramsey_nls)
  colnames(data_interest_ramsey_nls) =  c("output_gap", "inflation", "interest_rate","output_gap.l1","inflation.l1", "inflation.l2",
                                          "interest_rate.l1","interest_rate.l2","rer_gap", "rer_gap.l1", "inflation_ext.l2",
                                          "real_interest_gap","real_interest_gap.l1","real_interest_gap.l2","real_interest_gap.l3" ,"real_interest_gap.l4","real_interest_gap.l5", "real_interest_gap.l6",
                                          "us_output_gap","inflation_gap","brent_gap",
                                          "predicted_interest", "predicted_interest.2", "predicted_interest.3", "predicted_interest.4")
  
  reset_regression_interest = resettest(interest_rate ~ 1+inflation_gap+output_gap
                                          +interest_rate.l1 
                                          ,data = data_interest_ramsey_nls, type = "fitted", power = 1:4) 
  reset_regression_interest$p.value
  reset_regression_interest_nls = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+
                                             b4*x4+b5*x5+b6*x6+b7*x7
                                             -rho*b1*x8-rho*b2*x9-rho*b3*x10
                                          -rho*b4*x11-rho*b5*x12-rho*b6*x13-rho*b7*x14,
                                          start = list(rho = 0.03, b0 = 0.0, b1 = 0.9, b2 = 0.4, b3= 0.15, b4 = 0, b5 = 0, b6 = 0, b7 = 0),
                                          data = list(y = data_interest_ramsey_nls$interest_rate[2:obs_ramsey], y_l = data_interest_ramsey_nls$interest_rate[1:obs_ramsey-1], x1 = data_interest_ramsey_nls$interest_rate.l1[2:obs_ramsey],
                                                      x2 = data_interest_ramsey_nls$inflation_gap[2:obs_ramsey], x3 = data_interest_ramsey_nls$output_gap[2:obs_ramsey], 
                                                      x4 = data_interest_ramsey_nls$predicted_interest[2:obs_ramsey], x5 = data_interest_ramsey_nls$predicted_interest.2[2:obs_ramsey], x6 = data_interest_ramsey_nls$predicted_interest.3[2:obs_ramsey], x7 = data_interest_ramsey_nls$predicted_interest.4[2:obs_ramsey],  
                                                      x8 = data_interest_ramsey_nls$interest_rate.l1[1:obs_ramsey-1], x9 = data_interest_ramsey_nls$inflation_gap[1:obs_ramsey-1], x10 = data_interest_ramsey_nls$output_gap[1:obs_ramsey-1],
                                                      x11 = data_interest_ramsey_nls$predicted_interest[1:obs_ramsey-1],x12 = data_interest_ramsey_nls$predicted_interest.2[1:obs_ramsey-1], x13 = data_interest_ramsey_nls$predicted_interest.3[1:obs_ramsey-1], x14 = data_interest_ramsey_nls$predicted_interest.4[1:obs_ramsey-1])
                                          ,trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200)))
  
  summary(reset_regression_interest_nls)
  predicted_interest_nls_ramsey = predict(reset_regression_interest_nls)   
  ssr_interest_nls_ramsey <- sum((predicted_interest_nls_ramsey - data_interest_ramsey_nls$interest_rate[2:obs_ramsey]) ^ 2)  ## residual sum of squares
  #F-teset over with H0: predicted inf powers = 0
  F_interest_ramsey <-((ssr_1_interest_nls - ssr_interest_nls_ramsey)/(1+length(reset_regression_interest_nls$m$getPars())-length(interest_eq_nls$m$getPars())))/(ssr_interest_nls_ramsey/(nrow(data_interest_ramsey_nls[2:obs_ramsey,])-length(reset_regression_interest_nls$m$getPars())))
  pf(F_interest_ramsey,1+length(reset_regression_interest_nls$m$getPars())-length(interest_eq_nls$m$getPars()),nrow(data_interest_ramsey_nls[2:obs_ramsey,])-length(reset_regression_interest_nls$m$getPars()), lower.tail = FALSE)
  #Chi-Square version of the RESET test: 
  Chi_interest_ramsey <-nrow(data_interest_ramsey_nls[2:obs_ramsey,])*((ssr_1_interest_nls - ssr_interest_nls_ramsey)/(ssr_1_interest_nls))
  print(Chi_interest_ramsey)
  pchisq(Chi_interest_ramsey,1+length(reset_regression_interest_nls$m$getPars())-length(interest_eq_nls$m$getPars()), lower.tail = FALSE)
  
  #Test over RAMSEY Equation: 
  residuals_interest_ramsey = nlsResiduals(reset_regression_interest_nls)$resi1[,2]
  
  #Normality test
  jarque.bera.test(residuals(reset_regression_interest_nls)) #Normality   
  shapiro.test(residuals(reset_regression_interest_nls)) 
  #Heterocedasticity test (manual): 
  summary(lm(resid(reset_regression_interest_nls)^2 ~1+inflation_gap+output_gap
             +interest_rate.l1+predicted_interest.2+predicted_interest.3+predicted_interest.4,
             data = data_interest_ramsey_nls[2:obs_ramsey,])) #Homoskedasticiy at 5%
  # Serial Correlation: 
  lm_aux_bp = lm(resid(reset_regression_interest_nls)~ 1)
  bgtest(lm_aux_bp, order = 12)
  acf2(resid(reset_regression_interest_nls)) 
  
# ------------------ LM version test of Direct Non-linearity -------------------
# The test is conducted with the variables that belong to the explanatory variables in each case. Therefore, just the interactions with pot. state var are included.     
# LM is computed as: 1. Estimate model under linearity hypothesis and collect residuals and SSR (already done). 2. Estimate the model with the following interactions: 
# For each endogenous variable z_k: error = regular terms of the equation + sum of interactions of each regular term with z_k to the J power (unrestricted model)
# Compute F statistic. 
data_nls_lm = cbind(data_nls[2:obs,],residuals_inf_nls$resi1[,2], residuals_outputgap_nls$resi1[,2], residuals_interest_nls$resi1[,2])

names(data_nls_lm)

# Inflation Equation: 
data_nls_lm_inflation = as.data.frame(cbind(data_nls_lm$inflation, data_nls_lm$`residuals_inf_nls$resi1[, 2]`, data_nls_lm$inflation.l1, data_nls_lm$output_gap, data_nls_lm$output_gap.l1, data_nls_lm$real_interest_gap.l2,
                              data_nls_lm$rer_gap, data_nls_lm$inflation_ext.l2))
colnames(data_nls_lm_inflation) = c("inflation", "resid_inf", "inflation.l1", "output_gap", "output_gap.l1", "real_interest_gap.l2", "rer_gap", "inflation_ext.l2")
obs = nrow(data_nls_lm_inflation)
names_inf = colnames(data_nls_lm_inflation)
LM_inflation_direct = NULL
for(i in 3:ncol(data_nls_lm_inflation)){
   #From x1 to x12 are regular model variables (six normal + six lagged with rho interaction). 
   #From x13 to x18 are interaction terms with power = 1
   #From x19 to x24 are interaction terms with power = 2
   #From x25 to x30 are interaction terms with power = 3
   #From x31 to x48 are the interacted variables with rho, and lagged
   lm_inf_direct_inf_aux = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4
                           +b5*x5+b6*x6 
                           -rho*b1*x7-rho*b2*x8-rho*b3*x9-rho*b4*x10
                           -rho*b5*x11-rho*b6*x12
                           +b7*x13 + b8*x14 + b9*x15 + b10*x16 + b11*x17 + b12*x18
                           +b13*x19 + b14*x20 + b15*x21 + b16*x22 + b17*x23 + b18*x24
                           +b19*x25 + b20*x26 + b21*x27 + b22*x28 + b23*x29 + b24*x30
                           -rho*b7*x31 - rho*b8*x32 - rho*b9*x33 - rho*b10*x34 - rho*b11*x35 - rho*b12*x36
                           - rho*b13*x37 - rho*b14*x38 - rho*b15*x39 - rho*b16*x40 - rho*b17*x41 - rho*b18*x42
                           - rho*b19*x43 - rho*b20*x44 - rho*b21*x45 - rho*b22*x46 - rho*b23*x47 - rho*b24*x48
                           , 
                           start = list(rho = 0.01, b0 = 0.01, b1 = 0.67, b2 = -0.11, b3= 0.38, b4= -0.09, b5= 0.084, b6 = 0.02, 
                                        b7 = 0, b8 = 0, b9 = 0, b10 = 0, b11 = 0, b12 = 0, b13 = 0, b14 = 0, b15 = 0, b16 = 0, b17 = 0,
                                        b18 = 0, b19 = 0, b20 = 0, b21 = 0, b22 = 0, b23 = 0, b24 = 0),
                           data = list(y = data_nls_lm_inflation$resid_inf[2:obs], y_l = data_nls_lm_inflation$resid_inf[1:obs-1], 
                                       x1 = data_nls_lm_inflation$inflation.l1[2:obs], x2 = data_nls_lm_inflation$output_gap[2:obs], x3 = data_nls_lm_inflation$output_gap.l1[2:obs], 
                                       x4 = data_nls_lm_inflation$real_interest_gap.l2[2:obs], x5 = data_nls_lm_inflation$rer_gap[2:obs], x6 = data_nls_lm_inflation$inflation_ext.l2[2:obs],
                                       x7 = data_nls_lm_inflation$inflation.l1[1:obs-1], x8 = data_nls_lm_inflation$output_gap[1:obs-1], x9 = data_nls_lm_inflation$output_gap.l1[1:obs-1], 
                                       x10 = data_nls_lm_inflation$real_interest_gap.l2[1:obs-1], x11 = data_nls_lm_inflation$rer_gap[1:obs-1], x12 = data_nls_lm_inflation$inflation_ext.l2[1:obs-1],
                           x13 = data_nls_lm_inflation$inflation.l1[2:obs]*(data_nls_lm_inflation[2:obs,i]^1), x14 = data_nls_lm_inflation$output_gap[2:obs]*(data_nls_lm_inflation[2:obs,i]^1), x15 = data_nls_lm_inflation$output_gap.l1[2:obs]*(data_nls_lm_inflation[2:obs,i]^1), x16 = data_nls_lm_inflation$real_interest_gap.l2[2:obs]*(data_nls_lm_inflation[2:obs,i]^1), x17 = data_nls_lm_inflation$rer_gap[2:obs]*(data_nls_lm_inflation[2:obs,i]^1), x18 = data_nls_lm_inflation$inflation_ext.l2[2:obs]*(data_nls_lm_inflation[2:obs,i]^1),
                           x19 = data_nls_lm_inflation$inflation.l1[2:obs]*(data_nls_lm_inflation[2:obs,i]^2), x20 = data_nls_lm_inflation$output_gap[2:obs]*(data_nls_lm_inflation[2:obs,i]^2), x21 = data_nls_lm_inflation$output_gap.l1[2:obs]*(data_nls_lm_inflation[2:obs,i]^2), x22 = data_nls_lm_inflation$real_interest_gap.l2[2:obs]*(data_nls_lm_inflation[2:obs,i]^2), x23 = data_nls_lm_inflation$rer_gap[2:obs]*(data_nls_lm_inflation[2:obs,i]^2), x24 = data_nls_lm_inflation$inflation_ext.l2[2:obs]*(data_nls_lm_inflation[2:obs,i]^2),
                           x25 = data_nls_lm_inflation$inflation.l1[2:obs]*(data_nls_lm_inflation[2:obs,i]^3), x26 = data_nls_lm_inflation$output_gap[2:obs]*(data_nls_lm_inflation[2:obs,i]^3), x27 = data_nls_lm_inflation$output_gap.l1[2:obs]*(data_nls_lm_inflation[2:obs,i]^3), x28 = data_nls_lm_inflation$real_interest_gap.l2[2:obs]*(data_nls_lm_inflation[2:obs,i]^3), x29 = data_nls_lm_inflation$rer_gap[2:obs]*(data_nls_lm_inflation[2:obs,i]^3), x30 = data_nls_lm_inflation$inflation_ext.l2[2:obs]*(data_nls_lm_inflation[2:obs,i]^3),
                           x31 = data_nls_lm_inflation$inflation.l1[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^1), x32 = data_nls_lm_inflation$output_gap[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^1), x33 = data_nls_lm_inflation$output_gap.l1[1:obs-1]*(data_nls_lm_inflation[2:obs-1,i]^1), x34 = data_nls_lm_inflation$real_interest_gap.l2[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^1), x35 = data_nls_lm_inflation$rer_gap[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^1), x36 = data_nls_lm_inflation$inflation_ext.l2[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^1),
                           x37 = data_nls_lm_inflation$inflation.l1[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^2), x38 = data_nls_lm_inflation$output_gap[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^2), x39 = data_nls_lm_inflation$output_gap.l1[1:obs-1]*(data_nls_lm_inflation[2:obs-1,i]^2), x40 = data_nls_lm_inflation$real_interest_gap.l2[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^2), x41 = data_nls_lm_inflation$rer_gap[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^2), x42 = data_nls_lm_inflation$inflation_ext.l2[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^2),
                           x43 = data_nls_lm_inflation$inflation.l1[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^3), x44 = data_nls_lm_inflation$output_gap[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^3), x45 = data_nls_lm_inflation$output_gap.l1[1:obs-1]*(data_nls_lm_inflation[2:obs-1,i]^3), x46 = data_nls_lm_inflation$real_interest_gap.l2[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^3), x47 = data_nls_lm_inflation$rer_gap[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^3), x48 = data_nls_lm_inflation$inflation_ext.l2[1:obs-1]*(data_nls_lm_inflation[1:obs-1,i]^3)
                           
                           ), 
                           trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200), nls.control(warnOnly = TRUE)))
   
   predicted_lm_inf_direct = predict(lm_inf_direct_inf_aux)
   ssr_ur_lm = sum((predicted_lm_inf_direct - data_nls_lm_inflation$resid_inf[2:obs])^2)  ## residual sum of squares
   #Compute the F-test statistic
   F_test_statistic_lm = ((ssr_1_inf_nls - ssr_ur_lm)/(summary(inflation_eq_nls)$df[2]-summary(lm_inf_direct_inf_aux)$df[2]))/(ssr_ur_lm/(summary(lm_inf_direct_inf_aux)$df[2]))
   print(F_test_statistic_lm)  #Sugiere no linealidad al 5% respecto a la inflación rezagada 
   pval_F_lm = pf(F_test_statistic_lm,summary(inflation_eq_nls)$df[2]-summary(lm_inf_direct_inf_aux)$df[2],summary(lm_inf_direct_inf_aux)$df[2], lower.tail = FALSE)
   #Compute the Chi-test statistic
   Chi_test_statistic_lm = (obs-1)*((ssr_1_inf_nls - ssr_ur_lm)/ssr_1_inf_nls)
   print(Chi_test_statistic_lm)
   pval_chi_lm = pchisq(Chi_test_statistic_lm,summary(inflation_eq_nls)$df[2]-summary(lm_inf_direct_inf_aux)$df[2], lower.tail = FALSE)
   
   if(pval_F_lm <= 0.05 && pval_chi_lm <=0.05){test_lm = "Nonlinear"} else {test_lm = "Linear"}
   #Wald tests for selecting the STR functional form: 
   
   h03 = wald.test(b = as.vector(lm_inf_direct_inf_aux$m$getPars()),vcov(lm_inf_direct_inf_aux),Terms = c(21:26), df = summary(lm_inf_direct_inf_aux)$df[2]) 
   h02 = wald.test(b = as.vector(lm_inf_direct_inf_aux$m$getPars()),vcov(lm_inf_direct_inf_aux),Terms = c(15:19), df = summary(lm_inf_direct_inf_aux)$df[2]) 
   h01 = wald.test(b = as.vector(lm_inf_direct_inf_aux$m$getPars()),vcov(lm_inf_direct_inf_aux),Terms = c(9:14), df = summary(lm_inf_direct_inf_aux)$df[2]) 
      #H03 results:   
      F_wald_h03 = h03$result$Ftest[1] 
      pval_wald_h03_F = h03$result$Ftest[4]
      Chi_wald_h03 = h03$result$chi2[1] 
      pval_wald_h03_Chi = h03$result$chi2[3]
      #H02 results: 
      F_wald_h02 = h02$result$Ftest[1] 
      pval_wald_h02_F = h02$result$Ftest[4]
      Chi_wald_h02 = h02$result$chi2[1] 
      pval_wald_h02_Chi = h02$result$chi2[3]
      #H01 results:   
      F_wald_h01 = h01$result$Ftest[1] 
      pval_wald_h01_F = h01$result$Ftest[4]
      Chi_wald_h01 = h01$result$chi2[1] 
      pval_wald_h01_Chi = h01$result$chi2[3]
      if(pval_wald_h02_F <= pval_wald_h03_F && pval_wald_h02_F <= pval_wald_h01_F){str_model_F = "LSTR2-ESTR"} else{str_model_F = "LSTR1"}
      if(pval_wald_h02_Chi <= pval_wald_h03_Chi && pval_wald_h02_Chi <= pval_wald_h01_Chi){str_model_Chi = "LSTR2-ESTR"} else{str_model_Chi = "LSTR1"}
      final_aux = rbind(F_test_statistic_lm, pval_F_lm, Chi_test_statistic_lm, pval_chi_lm,test_lm,str_model_F,str_model_Chi)
      colnames(final_aux) = names_inf[i]
      rownames(final_aux) = c("F_test","F_test_pval","Chi_test", "Chi_test_pval","LM result", "STR model F", "STR model Chi")
      LM_inflation_direct = cbind(LM_inflation_direct, final_aux)
    }

write.csv(LM_inflation_direct, 'LM_inflation_direct.csv')

# Output Gap Equation: 
data_nls_lm_outputgap = as.data.frame(cbind(data_nls_lm$output_gap, data_nls_lm$`residuals_outputgap_nls$resi1[, 2]`, data_nls_lm$output_gap.l1, data_nls_lm$real_interest_gap.l2, data_nls_lm$us_output_gap, data_nls_lm$rer_gap,
                                            data_nls_lm$brent_gap))
colnames(data_nls_lm_outputgap) = c("output_gap", "resid_output_gap", "output_gap.l1","real_interest_gap.l2", "us_output_gap", "rer_gap", "brent_gap")
obs = nrow(data_nls_lm_outputgap)
names_output_gap = colnames(data_nls_lm_outputgap)
LM_output_gap_direct = NULL
for(i in 3:ncol(data_nls_lm_outputgap)){
   #From x1 to x10 are regular model variables (six normal + six lagged with rho interaction). 
   #From x11 to x15 are interaction terms with power = 1
   #From x16 to x20 are interaction terms with power = 2
   #From x21 to x25 are interaction terms with power = 3
   #From x26 to x48 are the interacted variables with rho, and lagged
   lm_output_direct_aux = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4+b5*x5
                                   -rho*b1*x6-rho*b2*x7-rho*b3*x8-rho*b4*x9-rho*b5*x10
                                  +b6*x11 + b7*x12 + b8*x13 + b9*x14 + b10*x15
                                  +b11*x16 + b12*x17 + b13*x18 + b14*x19 + b15*x20
                                  +b16*x21 + b17*x22 + b18*x23 + b19*x24 + b20*x25
                                  -rho*b6*x26 -rho*b7*x27 -rho*b8*x28 -rho*b9*x29 -rho*b10*x30
                                  -rho*b11*x31 -rho*b12*x32 -rho*b13*x33 -rho*b14*x34 -rho*b15*x35
                                  -rho*b16*x36 -rho*b17*x37 -rho*b18*x38 -rho*b19*x39 -rho*b20*x40
                                  ,
                                   start = list(rho = 0.03, b0 = -0.0015, b1 = 0.88, b2 = -0.023, b3= -0.084, b4= 0.23, b5 = 0.004, 
                                                b6 = 0, b7 = 0, b8 = 0, b9=0,b10=0,b11=0,b12=0,b13=0,b14=0,b15=0,b16=0,b17=0,b18=0,
                                                b19=0, b20=0),
                                   data = list(y = data_nls_lm_outputgap$resid_output_gap[2:obs], y_l = data_nls_lm_outputgap$resid_output_gap[1:obs-1], x1 = data_nls_lm_outputgap$output_gap.l1[2:obs],
                                               x2 = data_nls_lm_outputgap$real_interest_gap.l2[2:obs], x3 = data_nls_lm_outputgap$us_output_gap[2:obs], x4 = data_nls_lm_outputgap$rer_gap[2:obs], x5 = data_nls_lm_outputgap$brent_gap[2:obs], 
                                               x6 = data_nls_lm_outputgap$output_gap.l1[1:obs-1], x7 = data_nls_lm_outputgap$real_interest_gap.l2[1:obs-1], x8 = data_nls_lm_outputgap$us_output_gap[1:obs-1], x9 = data_nls_lm_outputgap$rer_gap[1:obs-1], x10 = data_nls_lm_outputgap$brent_gap[1:obs-1],
                                               x11 = data_nls_lm_outputgap$output_gap.l1[2:obs]*(data_nls_lm_outputgap[2:obs,i]^1), x12 = data_nls_lm_outputgap$real_interest_gap.l2[2:obs]*(data_nls_lm_outputgap[2:obs,i]^1), x13 = data_nls_lm_outputgap$us_output_gap[2:obs]*(data_nls_lm_outputgap[2:obs,i]^1), x14 = data_nls_lm_outputgap$rer_gap[2:obs]*(data_nls_lm_outputgap[2:obs,i]^1), x15 = data_nls_lm_outputgap$brent_gap[2:obs]*(data_nls_lm_outputgap[2:obs,i]^1),
                                               x16 = data_nls_lm_outputgap$output_gap.l1[2:obs]*(data_nls_lm_outputgap[2:obs,i]^2), x17 = data_nls_lm_outputgap$real_interest_gap.l2[2:obs]*(data_nls_lm_outputgap[2:obs,i]^2), x18 = data_nls_lm_outputgap$us_output_gap[2:obs]*(data_nls_lm_outputgap[2:obs,i]^2), x19 = data_nls_lm_outputgap$rer_gap[2:obs]*(data_nls_lm_outputgap[2:obs,i]^2), x20 = data_nls_lm_outputgap$brent_gap[2:obs]*(data_nls_lm_outputgap[2:obs,i]^2),
                                               x21 = data_nls_lm_outputgap$output_gap.l1[2:obs]*(data_nls_lm_outputgap[2:obs,i]^3), x22 = data_nls_lm_outputgap$real_interest_gap.l2[2:obs]*(data_nls_lm_outputgap[2:obs,i]^3), x23 = data_nls_lm_outputgap$us_output_gap[2:obs]*(data_nls_lm_outputgap[2:obs,i]^3), x24 = data_nls_lm_outputgap$rer_gap[2:obs]*(data_nls_lm_outputgap[2:obs,i]^3), x25 = data_nls_lm_outputgap$brent_gap[2:obs]*(data_nls_lm_outputgap[2:obs,i]^3),
                                               x26 = data_nls_lm_outputgap$output_gap.l1[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^1), x27 = data_nls_lm_outputgap$real_interest_gap.l2[1:obs-1]*(data_nls_lm_outputgap[2:obs-1,i]^1), x28 = data_nls_lm_outputgap$us_output_gap[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^1), x29 = data_nls_lm_outputgap$rer_gap[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^1), x30 = data_nls_lm_outputgap$brent_gap[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^1),
                                               x31 = data_nls_lm_outputgap$output_gap.l1[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^2), x32 = data_nls_lm_outputgap$real_interest_gap.l2[1:obs-1]*(data_nls_lm_outputgap[2:obs-1,i]^2), x33 = data_nls_lm_outputgap$us_output_gap[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^2), x34 = data_nls_lm_outputgap$rer_gap[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^2), x35 = data_nls_lm_outputgap$brent_gap[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^2),
                                               x36 = data_nls_lm_outputgap$output_gap.l1[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^3), x37 = data_nls_lm_outputgap$real_interest_gap.l2[1:obs-1]*(data_nls_lm_outputgap[2:obs-1,i]^3), x38 = data_nls_lm_outputgap$us_output_gap[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^3), x39 = data_nls_lm_outputgap$rer_gap[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^3), x40 = data_nls_lm_outputgap$brent_gap[1:obs-1]*(data_nls_lm_outputgap[1:obs-1,i]^3)
                                   )
                                   ,trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200), nls.control(warnOnly = TRUE)))
   
   predicted_lm_outputgap_direct = predict(lm_output_direct_aux)
   ssr_ur_lm = sum((predicted_lm_outputgap_direct - data_nls_lm_outputgap$resid_output_gap[2:obs])^2)  ## residual sum of squares
   #Compute the F-test statistic
   F_test_statistic_lm = ((ssr_1_outputgap_nls - ssr_ur_lm)/(summary(outputgap_eq_nls)$df[2]-summary(lm_output_direct_aux)$df[2]))/(ssr_ur_lm/(summary(lm_output_direct_aux)$df[2]))
   print(F_test_statistic_lm) 
   pval_F_lm = pf(F_test_statistic_lm,summary(outputgap_eq_nls)$df[2]-summary(lm_output_direct_aux)$df[2],summary(lm_output_direct_aux)$df[2], lower.tail = FALSE)
   #Compute the Chi-test statistic
   Chi_test_statistic_lm = (obs-1)*((ssr_1_outputgap_nls - ssr_ur_lm)/ssr_1_outputgap_nls)
   print(Chi_test_statistic_lm)
   pval_chi_lm = pchisq(Chi_test_statistic_lm,summary(outputgap_eq_nls)$df[2]-summary(lm_output_direct_aux)$df[2], lower.tail = FALSE)
   
   if(pval_F_lm <= 0.05 && pval_chi_lm <=0.05){test_lm = "Nonlinear"} else {test_lm = "Linear"}
   #Wald tests for selecting the STR functional form: 
   
   h03 = wald.test(b = as.vector(lm_output_direct_aux$m$getPars()),vcov(lm_output_direct_aux),Terms = c(18:22), df = summary(lm_output_direct_aux)$df[2]) 
   h02 = wald.test(b = as.vector(lm_output_direct_aux$m$getPars()),vcov(lm_output_direct_aux),Terms = c(13:17), df = summary(lm_output_direct_aux)$df[2]) 
   h01 = wald.test(b = as.vector(lm_output_direct_aux$m$getPars()),vcov(lm_output_direct_aux),Terms = c(8:12), df = summary(lm_output_direct_aux)$df[2]) 
   #H03 results:   
   F_wald_h03 = h03$result$Ftest[1] 
   pval_wald_h03_F = h03$result$Ftest[4]
   Chi_wald_h03 = h03$result$chi2[1] 
   pval_wald_h03_Chi = h03$result$chi2[3]
   #H02 results: 
   F_wald_h02 = h02$result$Ftest[1] 
   pval_wald_h02_F = h02$result$Ftest[4]
   Chi_wald_h02 = h02$result$chi2[1] 
   pval_wald_h02_Chi = h02$result$chi2[3]
   #H01 results:   
   F_wald_h01 = h01$result$Ftest[1] 
   pval_wald_h01_F = h01$result$Ftest[4]
   Chi_wald_h01 = h01$result$chi2[1] 
   pval_wald_h01_Chi = h01$result$chi2[3]
   if(pval_wald_h02_F <= pval_wald_h03_F && pval_wald_h02_F <= pval_wald_h01_F){str_model_F = "LSTR2-ESTR"} else{str_model_F = "LSTR1"}
   if(pval_wald_h02_Chi <= pval_wald_h03_Chi && pval_wald_h02_Chi <= pval_wald_h01_Chi){str_model_Chi = "LSTR2-ESTR"} else{str_model_Chi = "LSTR1"}
   final_aux = rbind(F_test_statistic_lm, pval_F_lm, Chi_test_statistic_lm, pval_chi_lm,test_lm,str_model_F,str_model_Chi)
   colnames(final_aux) = names_output_gap[i]
   rownames(final_aux) = c("F_test","F_test_pval","Chi_test", "Chi_test_pval","LM result", "STR model F", "STR model Chi")
   LM_output_gap_direct = cbind(LM_output_gap_direct, final_aux)
}

write.csv(LM_output_gap_direct, 'LM_outputgap_direct.csv')

# Interest rate Equation:  
data_nls_lm_interest = as.data.frame(cbind(data_nls_lm$interest_rate, data_nls_lm$`residuals_interest_nls$resi1[, 2]`, data_nls_lm$interest_rate.l1, data_nls_lm$inflation_gap, data_nls_lm$output_gap))
colnames(data_nls_lm_interest) = c("interest_rate", "resid_interest", "interest_rate.l1", "inflation_gap", "output_gap")
obs = nrow(data_nls_lm_interest)
names_interest = colnames(data_nls_lm_interest)
LM_interest_direct = NULL
for(i in 3:ncol(data_nls_lm_interest)){
   #From x1 to x6 are regular model variables (six normal + six lagged with rho interaction). 
   #From x7 to x9 are interaction terms with power = 1
   #From x10 to x12 are interaction terms with power = 2
   #From x13 to x15 are interaction terms with power = 3
   #From x16 to x24 are the interacted variables with rho, and lagged
   lm_interest_direct_aux = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+
                                     -rho*b1*x4-rho*b2*x5-rho*b3*x6
                                  +b4*x7 + b5*x8 + b6*x9
                                  +b7*x10 + b8*x11 + b9*x12
                                  +b10*x13 + b11*x14 + b12*x15
                                  -rho*b4*x16 -rho*b5*x17 -rho*b6*x18
                                  -rho*b7*x19 -rho*b8*x20 -rho*b9*x21
                                  -rho*b10*x22 -rho*b11*x23 -rho*b12*x24
                                  ,
                                  start = list(rho = 0.03, b0 = 0.0, b1 = 0.9, b2 = 0.4, b3= 0.15,
                                               b4 = 0,b5 = 0, b6 = 0, b7=0, b8 = 0, b9 = 0, b10=0, b11=0,b12=0),
                                  data = list(y = data_nls_lm_interest$resid_interest[2:obs], y_l = data_nls_lm_interest$resid_interest[1:obs-1], x1 = data_nls_lm_interest$interest_rate.l1[2:obs],
                                              x2 = data_nls_lm_interest$inflation_gap[2:obs], x3 = data_nls_lm_interest$output_gap[2:obs], x4 = data_nls_lm_interest$interest_rate.l1[1:obs-1], 
                                              x5 = data_nls_lm_interest$inflation_gap[1:obs-1], x6 = data_nls_lm_interest$output_gap[1:obs-1], 
                                              x7 = data_nls_lm_interest$interest_rate.l1[2:obs]*(data_nls_lm_interest[2:obs,i]^1), x8 = data_nls_lm_interest$inflation_gap[2:obs]*(data_nls_lm_interest[2:obs,i]^1), x9 = data_nls_lm_interest$output_gap[2:obs]*(data_nls_lm_interest[2:obs,i]^1),
                                              x10= data_nls_lm_interest$interest_rate.l1[2:obs]*(data_nls_lm_interest[2:obs,i]^2), x11 = data_nls_lm_interest$inflation_gap[2:obs]*(data_nls_lm_interest[2:obs,i]^2), x12 = data_nls_lm_interest$output_gap[2:obs]*(data_nls_lm_interest[2:obs,i]^2),
                                              x13= data_nls_lm_interest$interest_rate.l1[2:obs]*(data_nls_lm_interest[2:obs,i]^3), x14 = data_nls_lm_interest$inflation_gap[2:obs]*(data_nls_lm_interest[2:obs,i]^3), x15 = data_nls_lm_interest$output_gap[2:obs]*(data_nls_lm_interest[2:obs,i]^3),
                                              x16= data_nls_lm_interest$interest_rate.l1[1:obs-1]*(data_nls_lm_interest[1:obs-1,i]^1), x17 = data_nls_lm_interest$inflation_gap[1:obs-1]*(data_nls_lm_interest[1:obs-1,i]^1), x18 = data_nls_lm_interest$output_gap[1:obs-1]*(data_nls_lm_interest[1:obs-1,i]^1),
                                              x19= data_nls_lm_interest$interest_rate.l1[1:obs-1]*(data_nls_lm_interest[1:obs-1,i]^2), x20 = data_nls_lm_interest$inflation_gap[1:obs-1]*(data_nls_lm_interest[1:obs-1,i]^2), x21 = data_nls_lm_interest$output_gap[1:obs-1]*(data_nls_lm_interest[1:obs-1,i]^2),
                                              x22= data_nls_lm_interest$interest_rate.l1[1:obs-1]*(data_nls_lm_interest[1:obs-1,i]^3), x23 = data_nls_lm_interest$inflation_gap[1:obs-1]*(data_nls_lm_interest[1:obs-1,i]^3), x24 = data_nls_lm_interest$output_gap[1:obs-1]*(data_nls_lm_interest[1:obs-1,i]^3)
                                              
                                 )
                                 ,trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200), nls.control(warnOnly = TRUE)))
   
   predicted_lm_interest_direct = predict(lm_interest_direct_aux)
   ssr_ur_lm = sum((predicted_lm_interest_direct - data_nls_lm_interest$resid_interest[2:obs])^2)  ## residual sum of squares
   #Compute the F-test statistic
   F_test_statistic_lm = ((ssr_1_interest_nls - ssr_ur_lm)/(summary(interest_eq_nls)$df[2]-summary(lm_interest_direct_aux)$df[2]))/(ssr_ur_lm/(summary(lm_interest_direct_aux)$df[2]))
   print(F_test_statistic_lm) 
   pval_F_lm = pf(F_test_statistic_lm,summary(interest_eq_nls)$df[2]-summary(lm_interest_direct_aux)$df[2],summary(lm_interest_direct_aux)$df[2], lower.tail = FALSE)
   #Compute the Chi-test statistic
   Chi_test_statistic_lm = (obs-1)*((ssr_1_interest_nls - ssr_ur_lm)/ssr_1_interest_nls)
   print(Chi_test_statistic_lm)
   pval_chi_lm = pchisq(Chi_test_statistic_lm,summary(interest_eq_nls)$df[2]-summary(lm_interest_direct_aux)$df[2], lower.tail = FALSE)
   
   if(pval_F_lm <= 0.05 && pval_chi_lm <=0.05){test_lm = "Nonlinear"} else {test_lm = "Linear"}
   #Wald tests for selecting the STR functional form: 
   
   h03 = wald.test(b = as.vector(lm_interest_direct_aux$m$getPars()),vcov(lm_interest_direct_aux),Terms = c(12:14), df = summary(lm_interest_direct_aux)$df[2]) 
   h02 = wald.test(b = as.vector(lm_interest_direct_aux$m$getPars()),vcov(lm_interest_direct_aux),Terms = c(9:11), df = summary(lm_interest_direct_aux)$df[2]) 
   h01 = wald.test(b = as.vector(lm_interest_direct_aux$m$getPars()),vcov(lm_interest_direct_aux),Terms = c(6:8), df = summary(lm_interest_direct_aux)$df[2]) 
   #H03 results:   
   F_wald_h03 = h03$result$Ftest[1] 
   pval_wald_h03_F = h03$result$Ftest[4]
   Chi_wald_h03 = h03$result$chi2[1] 
   pval_wald_h03_Chi = h03$result$chi2[3]
   #H02 results: 
   F_wald_h02 = h02$result$Ftest[1] 
   pval_wald_h02_F = h02$result$Ftest[4]
   Chi_wald_h02 = h02$result$chi2[1] 
   pval_wald_h02_Chi = h02$result$chi2[3]
   #H01 results:   
   F_wald_h01 = h01$result$Ftest[1] 
   pval_wald_h01_F = h01$result$Ftest[4]
   Chi_wald_h01 = h01$result$chi2[1] 
   pval_wald_h01_Chi = h01$result$chi2[3]
   if(pval_wald_h02_F <= pval_wald_h03_F && pval_wald_h02_F <= pval_wald_h01_F){str_model_F = "LSTR2-ESTR"} else{str_model_F = "LSTR1"}
   if(pval_wald_h02_Chi <= pval_wald_h03_Chi && pval_wald_h02_Chi <= pval_wald_h01_Chi){str_model_Chi = "LSTR2-ESTR"} else{str_model_Chi = "LSTR1"}
   final_aux = rbind(F_test_statistic_lm, pval_F_lm, Chi_test_statistic_lm, pval_chi_lm,test_lm,str_model_F,str_model_Chi)
   colnames(final_aux) = names_interest[i]
   rownames(final_aux) = c("F_test","F_test_pval","Chi_test", "Chi_test_pval","LM result", "STR model F", "STR model Chi")
   LM_interest_direct = cbind(LM_interest_direct, final_aux)
}
write.csv(LM_interest_direct, 'LM_interest_direct.csv')

# ------------------ LM version test for State Dependency  ---------------------
# ------------------ and model selection (ESTAR, LSTAR1, LSTAR2)  ------------ #  
# LM is computed as: 1. Estimate model under linearity hypothesis and collect residuals. 2. Estimate the model with the following interactions: 
# For each endogenous variable z_k: error = regular terms of the equation +state variable to the 2-J power + sum of interactions of each regular term with z_k to the J power (unrestricted model)
# Compute F statistic. 
# Model selection is done using the Wald test approach
   
#Add state variables to data frame: 
   data_nls_lm = cbind(data_nls[2:nrow(data_nls),],residuals_inf_nls$resi1[,2], residuals_outputgap_nls$resi1[,2], residuals_interest_nls$resi1[,2])
   dates <- seq(as.Date("2002-12-01"), length=nrow(data_nls_lm), by="quarters")
   data_nls_lm = xts(data_nls_lm, order.by = dates)
   index(data_nls_lm)
   head(data_nls)
   colnames(data_nls_lm) = c("output_gap", "inflation", "interest_rate","output_gap.l1","inflation.l1", "inflation.l2",
                          "interest_rate.l1","interest_rate.l6","rer_gap", "rer_gap.l1", "inflation_ext.l2",
                          "real_interest_gap","real_interest_gap.l1","real_interest_gap.l2","real_interest_gap.l3","real_interest_gap.l4","real_interest_gap.l5","real_interest_gap.l6",
                          "us_output_gap", "inflation_gap", "brent_gap", "resid_inf", "resid_outputgap", "resid_interest")
   
   data_state = cbind.xts(data_xts$ICC,data_xts$IEC,data_xts$ICE, data_xts$IICV, data_xts$IICB, data_xts$IICVH, data_xts$ICCO, data_xts$ICCI,
                          data_xts$EPUC, data_xts$capacity_u_tc,data_xts$PIB, data_xts$PIB_d_index, data_xts$PIB_d_dep_index,
                          data_xts$NT_T, data_xts$investment_gdp, data_xts$consumption_y, data_xts$investment_yr, data_xts$agriculture, data_xts$minning, data_xts$manufacturing,
                          data_xts$retail, data_xts$construction, data_xts$services, data_xts$primary, data_xts$secondary, data_xts$terciary,data_xts$fiscal_budget_cop_yr, 
                          data_xts$gdp_pot_var, data_xts$TFP, data_xts$TFP_g, data_xts$rd_gdp, data_xts$hours_worked_g,data_xts$IRR,
                          data_xts$unempleyment_s,data_xts$tgp, data_xts$s_min, data_xts$s_comer,data_xts$informality, data_xts$labor_prod, 
                          data_xts$cash, data_xts$monetary_base, data_xts$m1, data_xts$m3, data_xts$CDT_pasivos, data_xts$Ahorros_pasivos, data_xts$Bonos_pasivos, data_xts$Bonos_pse, 
                          data_xts$Corriente_pasivos, data_xts$comercial_a, data_xts$comercial_icm, data_xts$comercial_icr, data_xts$comercial_cubr ,data_xts$consumo_a, data_xts$consumo_icm, data_xts$consumo_icr,data_xts$consumo_cubr,
                          data_xts$vivienda_a, data_xts$vivienda_icm, data_xts$vivienda_icr, data_xts$vivienda_cubr, data_xts$cartera_total_a, data_xts$vencida_total_a, data_xts$comercial_total, data_xts$consumo_total, data_xts$vivienda_total, 
                          data_xts$cubrimiento_total, data_xts$total_icm, data_xts$total_icr, data_xts$ICR_comer_desv, data_xts$ICR_cons_desv,data_xts$ICM_total_desv, data_xts$cartera_primaria_a, data_xts$cartera_secundaria_a, data_xts$cartera_terciaria_a, 
                          data_xts$cartera_primaria_p, data_xts$cartera_secundaria_p, data_xts$cartera_terciaria_p, data_xts$desem_c_plazo, data_xts$desem_comer_ord_plazo, data_xts$desem_comer_pf_plazo, data_xts$desem_total_plazo, 
                          data_xts$spread_plazo_consumo, data_xts$spread_plazo_ordinario, data_xts$spread_plazo_preferencial, data_xts$spread_plazo_total, data_xts$spread_ex_ante, data_xts$sprea_ex_post, data_xts$spread_e_prom_total_desv,
                          data_xts$Tasa_e_comer_desv, data_xts$Tasa_e_cons_desv, data_xts$Tasa_dep_e_desv, data_xts$NIM, data_xts$ROA, data_xts$caja_activos, data_xts$Cartera_activos, data_xts$Inversiones_activos,data_xts$Cartera_pse, data_xts$Cartera_PIB, 
                          data_xts$activos_patrimonio, data_xts$solvencia_total,data_xts$posicion_encaje, data_xts$M1_mult, data_xts$velocidad_m1, data_xts$IHH_Cartera, data_xts$IHH_Comercial, data_xts$IHH_Consumo, data_xts$IHH_Vivienda, data_xts$IHH_Ahorros, data_xts$IHH_CDT, 
                          data_xts$Lerner_promedio, data_xts$Lerner_desv, data_xts$deuda_ext_sr_lp,data_xts$deuda_ext_priv,data_xts$tasa_fija_deuda_ext, data_xts$interes_ex_post_deuda_ext,data_xts$COLCAP, data_xts$tes_1y, data_xts$tes_5y, data_xts$tes_10y, data_xts$`10yv1y`, 
                          data_xts$current_account_yr, data_xts$ied_yr, data_xts$ied_yr_primary, data_xts$ied_yr_secondary, data_xts$ied_yr_terciary,data_xts$fwd_deuda_ext, data_xts$trade_gdp, data_xts$world_gdp, data_xts$latam_gdp, data_xts$brent_real, data_xts$TRM,
                          data_xts$exp_inf_desv_pm4, data_xts$exp_trm_desv_pm4,data_xts$exp_dtf_desv_pm4, data_xts$exp_pib_desv_pm4, data_xts$exp_empleo_cp_pm4,data_xts$exp_liquidez_cp, data_xts$exp_credito_cp_pm4,   
                          data_xts$Voice_accuntability, data_xts$Pol_stability, data_xts$Gov_effective, data_xts$Regul_quality, data_xts$Rule_law, data_xts$Control_corruption, data_xts$governance_index, data_xts$Poverty_rate, data_xts$Gini,data_xts$offen_farc, data_xts$Homicidio_g, data_xts$Secuestro_g,data_xts$Hurto_g, 
                          lag(data_xts$output_gap_sm), lag(data_xts$inflation_sm), lag(data_xts$interest_rate_sm), lag(data_xts$output_gap_sv), lag(data_xts$inflation_sv), lag(data_xts$interest_rate_sv))  
   ncol(data_state)
   colnames(data_state) = c("ICC","IEC","ICE", "IICV", "IICB", "IICVH", "ICCO", "ICI", "EPUC", "capacity_u", "PIB", "PIB_d_index", "PIB_d_dep_index",
                            "NT_T", "investment_gdp", "consumption_yr", "investment_yr", "agriculture", "minning","manufacturing", "retail", "construction","services", "primary", "seconday", "terciary", "fiscal_budget", 
                            "gdp_pot", "TFP", "TFP_g", "rd_gdp", "hours_worked_g", "IRR", 
                            "unemployment", "tgp", "s_min", "s_comer","informality","labor_prod", 
                            "cash", "monetary_base","m1", "m3", "CDT_pasivos", "Ahorros_pasivos", "Bonos_pasivos", "Bonos_pse","Corriente_pasivos", "comercial_a", "comercial_icm", "comercial_icr","comcercial_cubr" ,"consumo_a","consumo_icm", "consumo_icr", "consumo_cubr",
                            "vivienda_a", "vivienda_icr", "vivienda_icm", "vivienda_cubr", "cartera_total_a", "vencida_total_a", "comercial_total", "consumo_total", "vivienda_total", 
                            "cubrimiento_total", "total_icm", "total_icr", "icr_comer_dessv", "icr_comer_desv", "icm_total_desv", "cartera_primaria_a", "cartera_secundaria_a", "cartera_terciaria_a", 
                            "cartera_primaria_p", "cartera_secundaria_p", "cartera_terciaria_p", "desem_c_plazo", "desem_comer_ord_plazo", "desem_comer_pf_plazo", "desem_total_plazo", 
                            "spread_plazo_consumo", "spread_plazo_ord", "spread_plazop_pf", "spread_plazo_total", "spread_ex_ante", "spread_ex_post", "spread_effectivo_desv",  
                            "tasa_e_comer_desv", "tasa_e_cons_desv", "tasa_e_dep_desv", "NIM", "ROA", "caja_activos", "Cartera_activos", "inversiones_activos", "cartera_pse", "cartera_PIB", 
                            "activos_patrimonio", "solvencia_total", "posicion_encaje", "m1_mult", "velocidad_mult",  "IHH_cartera", "IHH_comer", "IHH_cons", "IHH_vivienda", "IHH_ahorros", "IHH_CDT", 
                            "lerner_promedio", "lerner_desv", "deuda_ext_sr_lp", "deuda_ext_priv", "tasa_fija_deudaext", "interes_expost_deudaext", "COCLAP", "tes1y", "tes5y", "tes10y", "10YV1Y",
                            "current_account_yr", "ied_yr", "ied_yr_primary", "ied_yr_secndary", "ied_yr_terciary", "fwd_deuda_ext", "trad_gdp", "world_gdp", "latam_gdp", "brent_real", "TRM",
                            "exp_inf_desv", "exp_trm_desv", "exp_dtf_desv", "exp_pib_desv", "exp_empleo_cp", "exp_liquidez_cp", "exp_credito_cp", 
                            "voice_accountability", "pol_stability", "gov_effective", "regul_quality", "rule_law", "control_corruption", "governance_index", "poverty_rate", "Gini", "offen_farc", "Homicidio_g", "Secuestro_g", "Hurto_g", 
                            "output_gap_sm.l1", "inflation_sm.l1", "interest_sm.l1", "output_gap_sv.l1", "inflation_sv.l1", "interest_sv.l1")
   index(data_state)
   data_state_final = merge(data_nls_lm, data_state,join = "left") 
   index_dates_data_state_final = index(data_state_final)
   data_state_final = cbind.xts(data_state_final, xts(seq(1,nrow(data_state_final), 1),index_dates_data_state_final))
   names_data_state_final = c(colnames(data_nls_lm), colnames(data_state), "time_trend")
   colnames(data_state_final) = names_data_state_final
   colnames(data_state_final)
   
   #Orthogonalise potential state-dependency variables (except endogenous central tendency and dispersion measures)
   orthogonal_state = NULL
   for(i in eval(ncol(data_nls_lm)+1):eval(ncol(data_state_final)-7)){
      reg = lm(data_state_final[,i]~ data_state_final$output_gap +data_state_final$inflation + data_state_final$time_trend+1, data = data_state_final)
      aux = residuals(reg)
      orthogonal_state = cbind(orthogonal_state, aux)
   } 
   colnames(orthogonal_state) = c("ICC","IEC","ICE", "IICV", "IICB", "IICVH", "ICCO", "ICI", "EPUC", "capacity_u", "PIB", "PIB_d_index", "PIB_d_dep_index",
                                  "NT_T", "investment_gdp", "consumption_yr", "investment_yr", "agriculture", "minning","manufacturing", "retail", "construction","services", "primary", "seconday", "terciary", "fiscal_budget", 
                                  "gdp_pot", "TFP", "TFP_g", "rd_gdp", "hours_worked_g", "IRR", 
                                  "unemployment", "tgp", "s_min", "s_comer","informality","labor_prod", 
                                  "cash", "monetary_base","m1", "m3", "CDT_pasivos", "Ahorros_pasivos", "Bonos_pasivos", "Bonos_pse","Corriente_pasivos", "comercial_a", "comercial_icm", "comercial_icr","comcercial_cubr" ,"consumo_a","consumo_icm", "consumo_icr", "consumo_cubr",
                                  "vivienda_a", "vivienda_icr", "vivienda_icm", "vivienda_cubr", "cartera_total_a", "vencida_total_a", "comercial_total", "consumo_total", "vivienda_total", 
                                  "cubrimiento_total", "total_icm", "total_icr", "icr_comer_dessv", "icr_comer_desv", "icm_total_desv", "cartera_primaria_a", "cartera_secundaria_a", "cartera_terciaria_a", 
                                  "cartera_primaria_p", "cartera_secundaria_p", "cartera_terciaria_p", "desem_c_plazo", "desem_comer_ord_plazo", "desem_comer_pf_plazo", "desem_total_plazo", 
                                  "spread_plazo_consumo", "spread_plazo_ord", "spread_plazop_pf", "spread_plazo_total", "spread_ex_ante", "spread_ex_post", "spread_effectivo_desv",  
                                  "tasa_e_comer_desv", "tasa_e_cons_desv", "tasa_e_dep_desv", "NIM", "ROA", "caja_activos", "Cartera_activos", "inversiones_activos", "cartera_pse", "cartera_PIB", 
                                  "activos_patrimonio", "solvencia_total", "posicion_encaje", "m1_mult", "velocidad_mult",  "IHH_cartera", "IHH_comer", "IHH_cons", "IHH_vivienda", "IHH_ahorros", "IHH_CDT", 
                                  "lerner_promedio", "lerner_desv", "deuda_ext_sr_lp", "deuda_ext_priv", "tasa_fija_deudaext", "interes_expost_deudaext", "COCLAP", "tes1y", "tes5y", "tes10y", "10YV1Y",
                                  "current_account_yr", "ied_yr", "ied_yr_primary", "ied_yr_secndary", "ied_yr_terciary", "fwd_deuda_ext", "trad_gdp", "world_gdp", "latam_gdp", "brent_real", "TRM",
                                  "exp_inf_desv", "exp_trm_desv", "exp_dtf_desv", "exp_pib_desv", "exp_empleo_cp", "exp_liquidez_cp", "exp_credito_cp", 
                                  "voice_accountability", "pol_stability", "gov_effective", "regul_quality", "rule_law", "control_corruption", "governance_index", "poverty_rate", "Gini", "offen_farc", "Homicidio_g", "Secuestro_g", "Hurto_g" 
                                  )
# Inflation Equation: 
data_state_inflation = as.data.frame(cbind(data_state_final$inflation, data_state_final$resid_inf, data_state_final$inflation.l1, data_state_final$output_gap, data_state_final$output_gap.l1, data_state_final$real_interest_gap.l2,
                                               data_state_final$rer_gap, data_state_final$inflation_ext.l2))
data_state_inflation = cbind(data_state_inflation, orthogonal_state, data_state_final$output_gap_sm.l1, data_state_final$inflation_sm.l1,
                             data_state_final$interest_sm.l1, data_state_final$output_gap_sv.l1, data_state_final$inflation_sv.l1, data_state_final$interest_sv.l1)
obs = nrow(data_state_inflation)
names_inf = colnames(data_state_inflation)
LM_inflation_state = NULL   
for(i in 9:ncol(data_state_inflation)){
   #From x1 to x12 are regular model variables (six normal + six lagged with rho interaction). 
   #From x13 to x18 are interaction terms with power = 1
   #From x19 to x24 are interaction terms with power = 2
   #From x25 to x30 are interaction terms with power = 3
   #From x31 to x48 are the interacted variables with rho, and lagged
   #From x49 to x54 are state variables to the J power an its lagged interaction 
   lm_inf_direct_inf_aux = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4
                                   +b5*x5+b6*x6 
                                   -rho*b1*x7-rho*b2*x8-rho*b3*x9-rho*b4*x10
                                   -rho*b5*x11-rho*b6*x12
                                   +b7*x13 + b8*x14 + b9*x15 + b10*x16 + b11*x17 + b12*x18
                                   +b13*x19 + b14*x20 + b15*x21 + b16*x22 + b17*x23 + b18*x24
                                   +b19*x25 + b20*x26 + b21*x27 + b22*x28 + b23*x29 + b24*x30
                                   -rho*b7*x31 - rho*b8*x32 - rho*b9*x33 - rho*b10*x34 - rho*b11*x35 - rho*b12*x36
                                   - rho*b13*x37 - rho*b14*x38 - rho*b15*x39 - rho*b16*x40 - rho*b17*x41 - rho*b18*x42
                                   - rho*b19*x43 - rho*b20*x44 - rho*b21*x45 - rho*b22*x46 - rho*b23*x47 - rho*b24*x48
                                   +b25*x49 + b26*x50 + b27*x51
                                   -rho*b25*x52 -rho*b26*x53 -rho*b27*x54
                                   , 
                                   start = list(rho = 0.01, b0 = 0.01, b1 = 0.67, b2 = -0.11, b3= 0.38, b4= -0.09, b5= 0.084, b6 = 0.02, 
                                                b7 = 0, b8 = 0, b9 = 0, b10 = 0, b11 = 0, b12 = 0, b13 = 0, b14 = 0, b15 = 0, b16 = 0, b17 = 0,
                                                b18 = 0, b19 = 0, b20 = 0, b21 = 0, b22 = 0, b23 = 0, b24 = 0,
                                                b25 = 0, b26 = 0, b27 = 0),
                                   data = list(y = data_state_inflation$resid_inf[2:obs], y_l = data_state_inflation$resid_inf[1:obs-1], 
                                               x1 = data_state_inflation$inflation.l1[2:obs], x2 = data_state_inflation$output_gap[2:obs], x3 = data_state_inflation$output_gap.l1[2:obs], 
                                               x4 = data_state_inflation$real_interest_gap.l2[2:obs], x5 = data_state_inflation$rer_gap[2:obs], x6 = data_state_inflation$inflation_ext.l2[2:obs],
                                               x7 = data_state_inflation$inflation.l1[1:obs-1], x8 = data_state_inflation$output_gap[1:obs-1], x9 = data_state_inflation$output_gap.l1[1:obs-1], 
                                               x10 = data_state_inflation$real_interest_gap.l2[1:obs-1], x11 = data_state_inflation$rer_gap[1:obs-1], x12 = data_state_inflation$inflation_ext.l2[1:obs-1],
                                               x13 = data_state_inflation$inflation.l1[2:obs]*(data_state_inflation[2:obs,i]^1), x14 = data_state_inflation$output_gap[2:obs]*(data_state_inflation[2:obs,i]^1), x15 = data_state_inflation$output_gap.l1[2:obs]*(data_state_inflation[2:obs,i]^1), x16 = data_state_inflation$real_interest_gap.l2[2:obs]*(data_state_inflation[2:obs,i]^1), x17 = data_state_inflation$rer_gap[2:obs]*(data_state_inflation[2:obs,i]^1), x18 = data_state_inflation$inflation_ext.l2[2:obs]*(data_state_inflation[2:obs,i]^1),
                                               x19 = data_state_inflation$inflation.l1[2:obs]*(data_state_inflation[2:obs,i]^2), x20 = data_state_inflation$output_gap[2:obs]*(data_state_inflation[2:obs,i]^2), x21 = data_state_inflation$output_gap.l1[2:obs]*(data_state_inflation[2:obs,i]^2), x22 = data_state_inflation$real_interest_gap.l2[2:obs]*(data_state_inflation[2:obs,i]^2), x23 = data_state_inflation$rer_gap[2:obs]*(data_state_inflation[2:obs,i]^2), x24 = data_state_inflation$inflation_ext.l2[2:obs]*(data_state_inflation[2:obs,i]^2),
                                               x25 = data_state_inflation$inflation.l1[2:obs]*(data_state_inflation[2:obs,i]^3), x26 = data_state_inflation$output_gap[2:obs]*(data_state_inflation[2:obs,i]^3), x27 = data_state_inflation$output_gap.l1[2:obs]*(data_state_inflation[2:obs,i]^3), x28 = data_state_inflation$real_interest_gap.l2[2:obs]*(data_state_inflation[2:obs,i]^3), x29 = data_state_inflation$rer_gap[2:obs]*(data_state_inflation[2:obs,i]^3), x30 = data_state_inflation$inflation_ext.l2[2:obs]*(data_state_inflation[2:obs,i]^3),
                                               x31 = data_state_inflation$inflation.l1[1:obs-1]*(data_state_inflation[1:obs-1,i]^1), x32 = data_state_inflation$output_gap[1:obs-1]*(data_state_inflation[1:obs-1,i]^1), x33 = data_state_inflation$output_gap.l1[1:obs-1]*(data_state_inflation[2:obs-1,i]^1), x34 = data_state_inflation$real_interest_gap.l2[1:obs-1]*(data_state_inflation[1:obs-1,i]^1), x35 = data_state_inflation$rer_gap[1:obs-1]*(data_state_inflation[1:obs-1,i]^1), x36 = data_state_inflation$inflation_ext.l2[1:obs-1]*(data_state_inflation[1:obs-1,i]^1),
                                               x37 = data_state_inflation$inflation.l1[1:obs-1]*(data_state_inflation[1:obs-1,i]^2), x38 = data_state_inflation$output_gap[1:obs-1]*(data_state_inflation[1:obs-1,i]^2), x39 = data_state_inflation$output_gap.l1[1:obs-1]*(data_state_inflation[2:obs-1,i]^2), x40 = data_state_inflation$real_interest_gap.l2[1:obs-1]*(data_state_inflation[1:obs-1,i]^2), x41 = data_state_inflation$rer_gap[1:obs-1]*(data_state_inflation[1:obs-1,i]^2), x42 = data_state_inflation$inflation_ext.l2[1:obs-1]*(data_state_inflation[1:obs-1,i]^2),
                                               x43 = data_state_inflation$inflation.l1[1:obs-1]*(data_state_inflation[1:obs-1,i]^3), x44 = data_state_inflation$output_gap[1:obs-1]*(data_state_inflation[1:obs-1,i]^3), x45 = data_state_inflation$output_gap.l1[1:obs-1]*(data_state_inflation[2:obs-1,i]^3), x46 = data_state_inflation$real_interest_gap.l2[1:obs-1]*(data_state_inflation[1:obs-1,i]^3), x47 = data_state_inflation$rer_gap[1:obs-1]*(data_state_inflation[1:obs-1,i]^3), x48 = data_state_inflation$inflation_ext.l2[1:obs-1]*(data_state_inflation[1:obs-1,i]^3),
                                               x49 = data_state_inflation[2:obs,i]^1,x50 = data_state_inflation[2:obs,i]^2, x51 =data_state_inflation[2:obs,i]^3, 
                                               x52 = data_state_inflation[1:obs-1,i]^1,x53 = data_state_inflation[1:obs-1,i]^2, x54 =data_state_inflation[1:obs-1,i]^3
                                               ), 
                                   trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200), nls.control(warnOnly = TRUE)), silent= TRUE)
   if(try(summary(lm_inf_direct_inf_aux), silent = TRUE)[2] %in% 'try-error' || class(lm_inf_direct_inf_aux) %in% 'try-error'){next} else {
   #summary(lm_inf_direct_inf_aux)
   predicted_lm_inf_direct = predict(lm_inf_direct_inf_aux)
   ssr_ur_lm = sum((predicted_lm_inf_direct - data_nls_lm_inflation$resid_inf[2:obs])^2)  ## residual sum of squares
   #Compute the F-test statistic
   F_test_statistic_lm = ((ssr_1_inf_nls - ssr_ur_lm)/(summary(inflation_eq_nls)$df[2]-summary(lm_inf_direct_inf_aux)$df[2]))/(ssr_ur_lm/(summary(lm_inf_direct_inf_aux)$df[2]))
   print(F_test_statistic_lm)  #Sugiere no linealidad al 5% respecto a la inflación rezagada 
   pval_F_lm = pf(F_test_statistic_lm,summary(inflation_eq_nls)$df[2]-summary(lm_inf_direct_inf_aux)$df[2],summary(lm_inf_direct_inf_aux)$df[2], lower.tail = FALSE)
   #Compute the Chi-test statistic
   Chi_test_statistic_lm = (obs-1)*((ssr_1_inf_nls - ssr_ur_lm)/ssr_1_inf_nls)
   print(Chi_test_statistic_lm)
   pval_chi_lm = pchisq(Chi_test_statistic_lm,summary(inflation_eq_nls)$df[2]-summary(lm_inf_direct_inf_aux)$df[2], lower.tail = FALSE)
   
   if(pval_F_lm <= 0.05 && pval_chi_lm <=0.05){test_lm = "Nonlinear"} else {test_lm = "Linear"}
   #Wald tests for selecting the STR functional form: 
   
   h03 = wald.test(b = as.vector(lm_inf_direct_inf_aux$m$getPars()),vcov(lm_inf_direct_inf_aux),Terms = c(29,21:26), df = summary(lm_inf_direct_inf_aux)$df[2]) 
   h02 = wald.test(b = as.vector(lm_inf_direct_inf_aux$m$getPars()),vcov(lm_inf_direct_inf_aux),Terms = c(28,15:19), df = summary(lm_inf_direct_inf_aux)$df[2]) 
   h01 = wald.test(b = as.vector(lm_inf_direct_inf_aux$m$getPars()),vcov(lm_inf_direct_inf_aux),Terms = c(27,9:14), df = summary(lm_inf_direct_inf_aux)$df[2]) 
   #H03 results:   
   F_wald_h03 = h03$result$Ftest[1] 
   pval_wald_h03_F = h03$result$Ftest[4]
   Chi_wald_h03 = h03$result$chi2[1] 
   pval_wald_h03_Chi = h03$result$chi2[3]
   #H02 results: 
   F_wald_h02 = h02$result$Ftest[1] 
   pval_wald_h02_F = h02$result$Ftest[4]
   Chi_wald_h02 = h02$result$chi2[1] 
   pval_wald_h02_Chi = h02$result$chi2[3]
   #H01 results:   
   F_wald_h01 = h01$result$Ftest[1] 
   pval_wald_h01_F = h01$result$Ftest[4]
   Chi_wald_h01 = h01$result$chi2[1] 
   pval_wald_h01_Chi = h01$result$chi2[3]
   if(pval_wald_h02_F <= pval_wald_h03_F && pval_wald_h02_F <= pval_wald_h01_F){str_model_F = "LSTR2-ESTR"} else{str_model_F = "LSTR1"}
   if(pval_wald_h02_Chi <= pval_wald_h03_Chi && pval_wald_h02_Chi <= pval_wald_h01_Chi){str_model_Chi = "LSTR2-ESTR"} else{str_model_Chi = "LSTR1"}
   final_aux = rbind(F_test_statistic_lm, pval_F_lm, Chi_test_statistic_lm, pval_chi_lm,test_lm,str_model_F,str_model_Chi)
   colnames(final_aux) = names_inf[i]
   rownames(final_aux) = c("F_test","F_test_pval","Chi_test", "Chi_test_pval","LM result", "STR model F", "STR model Chi")
   LM_inflation_state = cbind(LM_inflation_state, final_aux)}
}  
write.csv(LM_inflation_state, 'LM_inflation_state.csv')
   
# Output Gap Equation: 
data_state_outputgap = as.data.frame(cbind(data_state_final$output_gap, data_state_final$resid_outputgap, data_state_final$output_gap.l1, data_state_final$real_interest_gap.l2, data_state_final$us_output_gap, data_state_final$rer_gap,
                                           data_state_final$brent_gap))
data_state_outputgap = cbind(data_state_outputgap, orthogonal_state, data_state_final$output_gap_sm.l1, data_state_final$inflation_sm.l1,
                             data_state_final$interest_sm.l1, data_state_final$output_gap_sv.l1, data_state_final$inflation_sv.l1, data_state_final$interest_sv.l1)
obs = nrow(data_state_outputgap)
names_outputgap = colnames(data_state_outputgap)
LM_outputgap_state = NULL   

for(i in 8:ncol(data_state_outputgap)){
   #From x1 to x10 are regular model variables (six normal + six lagged with rho interaction). 
   #From x11 to x15 are interaction terms with power = 1
   #From x16 to x20 are interaction terms with power = 2
   #From x21 to x25 are interaction terms with power = 3
   #From x26 to x40 are the interacted variables with rho, and lagged
   #From x41 to x46 are the state variable to the J power + lagged with rho interaction.  
   lm_output_direct_aux = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4+b5*x5
                                  -rho*b1*x6-rho*b2*x7-rho*b3*x8-rho*b4*x9-rho*b5*x10
                                  +b6*x11 + b7*x12 + b8*x13 + b9*x14 + b10*x15
                                  +b11*x16 + b12*x17 + b13*x18 + b14*x19 + b15*x20
                                  +b16*x21 + b17*x22 + b18*x23 + b19*x24 + b20*x25
                                  -rho*b6*x26 -rho*b7*x27 -rho*b8*x28 -rho*b9*x29 -rho*b10*x30
                                  -rho*b11*x31 -rho*b12*x32 -rho*b13*x33 -rho*b14*x34 -rho*b15*x35
                                  -rho*b16*x36 -rho*b17*x37 -rho*b18*x38 -rho*b19*x39 -rho*b20*x40
                                  +b21*x41 + b22*x42 + b23*x43
                                  -rho*b21*x44 - rho*b22*x45 -rho*b23*x46
                                  ,
                                  start = list(rho = 0.03, b0 = -0.0015, b1 = 0.88, b2 = -0.023, b3= -0.084, b4= 0.23, b5 = 0.004, 
                                               b6 = 0, b7 = 0, b8 = 0, b9=0,b10=0,b11=0,b12=0,b13=0,b14=0,b15=0,b16=0,b17=0,b18=0,
                                               b19=0, b20=0, b21 = 0, b22 = 0, b23 = 0),
                                  data = list(y = data_state_outputgap$resid_outputgap[2:obs], y_l = data_state_outputgap$resid_outputgap[1:obs-1], x1 = data_state_outputgap$output_gap.l1[2:obs],
                                              x2 = data_state_outputgap$real_interest_gap.l2[2:obs], x3 = data_state_outputgap$us_output_gap[2:obs], x4 = data_state_outputgap$rer_gap[2:obs], x5 = data_state_outputgap$brent_gap[2:obs], 
                                              x6 = data_state_outputgap$output_gap.l1[1:obs-1], x7 = data_state_outputgap$real_interest_gap.l2[1:obs-1], x8 = data_state_outputgap$us_output_gap[1:obs-1], x9 = data_state_outputgap$rer_gap[1:obs-1], x10 = data_state_outputgap$brent_gap[1:obs-1],
                                              x11 = data_state_outputgap$output_gap.l1[2:obs]*(data_state_outputgap[2:obs,i]^1), x12 = data_state_outputgap$real_interest_gap.l2[2:obs]*(data_state_outputgap[2:obs,i]^1), x13 = data_state_outputgap$us_output_gap[2:obs]*(data_state_outputgap[2:obs,i]^1), x14 = data_state_outputgap$rer_gap[2:obs]*(data_state_outputgap[2:obs,i]^1), x15 = data_state_outputgap$brent_gap[2:obs]*(data_state_outputgap[2:obs,i]^1),
                                              x16 = data_state_outputgap$output_gap.l1[2:obs]*(data_state_outputgap[2:obs,i]^2), x17 = data_state_outputgap$real_interest_gap.l2[2:obs]*(data_state_outputgap[2:obs,i]^2), x18 = data_state_outputgap$us_output_gap[2:obs]*(data_state_outputgap[2:obs,i]^2), x19 = data_state_outputgap$rer_gap[2:obs]*(data_state_outputgap[2:obs,i]^2), x20 = data_state_outputgap$brent_gap[2:obs]*(data_state_outputgap[2:obs,i]^2),
                                              x21 = data_state_outputgap$output_gap.l1[2:obs]*(data_state_outputgap[2:obs,i]^3), x22 = data_state_outputgap$real_interest_gap.l2[2:obs]*(data_state_outputgap[2:obs,i]^3), x23 = data_state_outputgap$us_output_gap[2:obs]*(data_state_outputgap[2:obs,i]^3), x24 = data_state_outputgap$rer_gap[2:obs]*(data_state_outputgap[2:obs,i]^3), x25 = data_state_outputgap$brent_gap[2:obs]*(data_state_outputgap[2:obs,i]^3),
                                              x26 = data_state_outputgap$output_gap.l1[1:obs-1]*(data_state_outputgap[1:obs-1,i]^1), x27 = data_state_outputgap$real_interest_gap.l2[1:obs-1]*(data_state_outputgap[2:obs-1,i]^1), x28 = data_state_outputgap$us_output_gap[1:obs-1]*(data_state_outputgap[1:obs-1,i]^1), x29 = data_state_outputgap$rer_gap[1:obs-1]*(data_state_outputgap[1:obs-1,i]^1), x30 = data_state_outputgap$brent_gap[1:obs-1]*(data_state_outputgap[1:obs-1,i]^1),
                                              x31 = data_state_outputgap$output_gap.l1[1:obs-1]*(data_state_outputgap[1:obs-1,i]^2), x32 = data_state_outputgap$real_interest_gap.l2[1:obs-1]*(data_state_outputgap[2:obs-1,i]^2), x33 = data_state_outputgap$us_output_gap[1:obs-1]*(data_state_outputgap[1:obs-1,i]^2), x34 = data_state_outputgap$rer_gap[1:obs-1]*(data_state_outputgap[1:obs-1,i]^2), x35 = data_state_outputgap$brent_gap[1:obs-1]*(data_state_outputgap[1:obs-1,i]^2),
                                              x36 = data_state_outputgap$output_gap.l1[1:obs-1]*(data_state_outputgap[1:obs-1,i]^3), x37 = data_state_outputgap$real_interest_gap.l2[1:obs-1]*(data_state_outputgap[2:obs-1,i]^3), x38 = data_state_outputgap$us_output_gap[1:obs-1]*(data_state_outputgap[1:obs-1,i]^3), x39 = data_state_outputgap$rer_gap[1:obs-1]*(data_state_outputgap[1:obs-1,i]^3), x40 = data_state_outputgap$brent_gap[1:obs-1]*(data_state_outputgap[1:obs-1,i]^3),
                                              x41 = data_state_outputgap[2:obs,i]^1, x42 = data_state_outputgap[2:obs,i]^2, x43 = data_state_outputgap[2:obs,i]^3, 
                                              x44 = data_state_outputgap[1:obs-1,i]^1, x45 = data_state_outputgap[1:obs-1,i]^2, x46 = data_state_outputgap[1:obs-1,i]^3
                                  )
                                  ,trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200), nls.control(warnOnly = TRUE)), silent = TRUE)
   if(try(summary(lm_output_direct_aux), silent = TRUE)[2] %in% 'try-error' || class(lm_output_direct_aux) %in% 'try-error'){next} else {
   
   predicted_lm_outputgap_direct = predict(lm_output_direct_aux)
   ssr_ur_lm = sum((predicted_lm_outputgap_direct - data_nls_lm_outputgap$resid_output_gap[2:obs])^2)  ## residual sum of squares
   #Compute the F-test statistic
   F_test_statistic_lm = ((ssr_1_outputgap_nls - ssr_ur_lm)/(summary(outputgap_eq_nls)$df[2]-summary(lm_output_direct_aux)$df[2]))/(ssr_ur_lm/(summary(lm_output_direct_aux)$df[2]))
   print(F_test_statistic_lm) 
   pval_F_lm = pf(F_test_statistic_lm,summary(outputgap_eq_nls)$df[2]-summary(lm_output_direct_aux)$df[2],summary(lm_output_direct_aux)$df[2], lower.tail = FALSE)
   #Compute the Chi-test statistic
   Chi_test_statistic_lm = (obs-1)*((ssr_1_outputgap_nls - ssr_ur_lm)/ssr_1_outputgap_nls)
   print(Chi_test_statistic_lm)
   pval_chi_lm = pchisq(Chi_test_statistic_lm,summary(outputgap_eq_nls)$df[2]-summary(lm_output_direct_aux)$df[2], lower.tail = FALSE)
   
   if(pval_F_lm <= 0.05 && pval_chi_lm <=0.05){test_lm = "Nonlinear"} else {test_lm = "Linear"}
   #Wald tests for selecting the STR functional form: 
   
   h03 = wald.test(b = as.vector(lm_output_direct_aux$m$getPars()),vcov(lm_output_direct_aux),Terms = c(25,18:22), df = summary(lm_output_direct_aux)$df[2]) 
   h02 = wald.test(b = as.vector(lm_output_direct_aux$m$getPars()),vcov(lm_output_direct_aux),Terms = c(24,13:17), df = summary(lm_output_direct_aux)$df[2]) 
   h01 = wald.test(b = as.vector(lm_output_direct_aux$m$getPars()),vcov(lm_output_direct_aux),Terms = c(23,8:12), df = summary(lm_output_direct_aux)$df[2]) 
   #H03 results:   
   F_wald_h03 = h03$result$Ftest[1] 
   pval_wald_h03_F = h03$result$Ftest[4]
   Chi_wald_h03 = h03$result$chi2[1] 
   pval_wald_h03_Chi = h03$result$chi2[3]
   #H02 results: 
   F_wald_h02 = h02$result$Ftest[1] 
   pval_wald_h02_F = h02$result$Ftest[4]
   Chi_wald_h02 = h02$result$chi2[1] 
   pval_wald_h02_Chi = h02$result$chi2[3]
   #H01 results:   
   F_wald_h01 = h01$result$Ftest[1] 
   pval_wald_h01_F = h01$result$Ftest[4]
   Chi_wald_h01 = h01$result$chi2[1] 
   pval_wald_h01_Chi = h01$result$chi2[3]
   if(pval_wald_h02_F <= pval_wald_h03_F && pval_wald_h02_F <= pval_wald_h01_F){str_model_F = "LSTR2-ESTR"} else{str_model_F = "LSTR1"}
   if(pval_wald_h02_Chi <= pval_wald_h03_Chi && pval_wald_h02_Chi <= pval_wald_h01_Chi){str_model_Chi = "LSTR2-ESTR"} else{str_model_Chi = "LSTR1"}
   final_aux = rbind(F_test_statistic_lm, pval_F_lm, Chi_test_statistic_lm, pval_chi_lm,test_lm,str_model_F,str_model_Chi)
   colnames(final_aux) = names_outputgap[i]
   rownames(final_aux) = c("F_test","F_test_pval","Chi_test", "Chi_test_pval","LM result", "STR model F", "STR model Chi")
   LM_outputgap_state = cbind(LM_outputgap_state, final_aux)}
}
write.csv(LM_outputgap_state, 'LM_outputgap_state.csv')


# Interest rate Equation:    
data_state_interest = as.data.frame(cbind(data_state_final$interest_rate, data_state_final$resid_interest, data_state_final$interest_rate.l1, data_state_final$inflation_gap, data_state_final$output_gap))
data_state_interest = cbind(data_state_interest, orthogonal_state, data_state_final$output_gap_sm.l1, data_state_final$inflation_sm.l1,
                             data_state_final$interest_sm.l1, data_state_final$output_gap_sv.l1, data_state_final$inflation_sv.l1, data_state_final$interest_sv.l1)
obs = nrow(data_state_interest)
names_interest = colnames(data_state_interest)
LM_interest_state = NULL   
for(i in 6:ncol(data_state_interest)){
   #From x1 to x6 are regular model variables (six normal + six lagged with rho interaction). 
   #From x7 to x9 are interaction terms with power = 1
   #From x10 to x12 are interaction terms with power = 2
   #From x13 to x15 are interaction terms with power = 3
   #From x16 to x24 are the interacted variables with rho, and lagged
   #From x25 to x30 are the state variables to the J power and the lagged with rho interaction
   lm_interest_direct_aux = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+
                                        -rho*b1*x4-rho*b2*x5-rho*b3*x6
                                      +b4*x7 + b5*x8 + b6*x9
                                      +b7*x10 + b8*x11 + b9*x12
                                      +b10*x13 + b11*x14 + b12*x15
                                      -rho*b4*x16 -rho*b5*x17 -rho*b6*x18
                                      -rho*b7*x19 -rho*b8*x20 -rho*b9*x21
                                      -rho*b10*x22 -rho*b11*x23 -rho*b12*x24
                                      +b13*x25 + b14*x26 + b15*x27
                                      -rho*b13*x28 - rho*b14*x29 - rho*b15*x30
                                      ,
                                      start = list(rho = 0.03, b0 = 0.0, b1 = 0.9, b2 = 0.4, b3= 0.15,
                                                   b4 = 0,b5 = 0, b6 = 0, b7=0, b8 = 0, b9 = 0, b10=0, b11=0,b12=0,
                                                   b13 = 0,b14 = 0, b15 = 0),
                                      data = list(y = data_state_interest$resid_interest[2:obs], y_l = data_state_interest$resid_interest[1:obs-1], x1 = data_state_interest$interest_rate.l1[2:obs],
                                                  x2 = data_state_interest$inflation_gap[2:obs], x3 = data_state_interest$output_gap[2:obs], x4 = data_state_interest$interest_rate.l1[1:obs-1], 
                                                  x5 = data_state_interest$inflation_gap[1:obs-1], x6 = data_state_interest$output_gap[1:obs-1], 
                                                  x7 = data_state_interest$interest_rate.l1[2:obs]*(data_state_interest[2:obs,i]^1), x8 = data_state_interest$inflation_gap[2:obs]*(data_state_interest[2:obs,i]^1), x9 = data_state_interest$output_gap[2:obs]*(data_state_interest[2:obs,i]^1),
                                                  x10= data_state_interest$interest_rate.l1[2:obs]*(data_state_interest[2:obs,i]^2), x11 = data_state_interest$inflation_gap[2:obs]*(data_state_interest[2:obs,i]^2), x12 = data_state_interest$output_gap[2:obs]*(data_state_interest[2:obs,i]^2),
                                                  x13= data_state_interest$interest_rate.l1[2:obs]*(data_state_interest[2:obs,i]^3), x14 = data_state_interest$inflation_gap[2:obs]*(data_state_interest[2:obs,i]^3), x15 = data_state_interest$output_gap[2:obs]*(data_state_interest[2:obs,i]^3),
                                                  x16= data_state_interest$interest_rate.l1[1:obs-1]*(data_state_interest[1:obs-1,i]^1), x17 = data_state_interest$inflation_gap[1:obs-1]*(data_state_interest[1:obs-1,i]^1), x18 = data_state_interest$output_gap[1:obs-1]*(data_state_interest[1:obs-1,i]^1),
                                                  x19= data_state_interest$interest_rate.l1[1:obs-1]*(data_state_interest[1:obs-1,i]^2), x20 = data_state_interest$inflation_gap[1:obs-1]*(data_state_interest[1:obs-1,i]^2), x21 = data_state_interest$output_gap[1:obs-1]*(data_state_interest[1:obs-1,i]^2),
                                                  x22= data_state_interest$interest_rate.l1[1:obs-1]*(data_state_interest[1:obs-1,i]^3), x23 = data_state_interest$inflation_gap[1:obs-1]*(data_state_interest[1:obs-1,i]^3), x24 = data_state_interest$output_gap[1:obs-1]*(data_state_interest[1:obs-1,i]^3),
                                                  x25 = data_state_interest[2:obs,i]^1, x26 = data_state_interest[2:obs,i]^2, x27 = data_state_interest[2:obs,i]^3, 
                                                  x28 = data_state_interest[1:obs-1,i]^1, x29 = data_state_interest[1:obs,i-1]^2, x30 = data_state_interest[1:obs-1,i]^3 
                                                  
                                      )
                                      ,trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200), nls.control(warnOnly = TRUE)), silent = TRUE) 
   if(try(summary(lm_interest_direct_aux), silent = TRUE)[2] %in% 'try-error' || class(lm_interest_direct_aux) %in% 'try-error') {next} else 
   {predicted_lm_interest_direct = predict(lm_interest_direct_aux)
   ssr_ur_lm = sum((predicted_lm_interest_direct - data_state_interest$resid_interest[2:obs])^2)  ## residual sum of squares
   #Compute the F-test statistic
   F_test_statistic_lm = ((ssr_1_interest_nls - ssr_ur_lm)/(summary(interest_eq_nls)$df[2]-summary(lm_interest_direct_aux)$df[2]))/(ssr_ur_lm/(summary(lm_interest_direct_aux)$df[2]))
   print(F_test_statistic_lm) 
   pval_F_lm = pf(F_test_statistic_lm,summary(interest_eq_nls)$df[2]-summary(lm_interest_direct_aux)$df[2],summary(lm_interest_direct_aux)$df[2], lower.tail = FALSE)
   #Compute the Chi-test statistic
   Chi_test_statistic_lm = (obs-1)*((ssr_1_interest_nls - ssr_ur_lm)/ssr_1_interest_nls)
   print(Chi_test_statistic_lm)
   pval_chi_lm = pchisq(Chi_test_statistic_lm,summary(interest_eq_nls)$df[2]-summary(lm_interest_direct_aux)$df[2], lower.tail = FALSE)
   
   if(pval_F_lm <= 0.05 && pval_chi_lm <=0.05){test_lm = "Nonlinear"} else {test_lm = "Linear"}
   #Wald tests for selecting the STR functional form: 

   h03 = try(wald.test(b = as.vector(lm_interest_direct_aux$m$getPars()),vcov(lm_interest_direct_aux),Terms = c(17, 12:14), df = summary(lm_interest_direct_aux)$df[2]), silent = TRUE) 
   h02 = try(wald.test(b = as.vector(lm_interest_direct_aux$m$getPars()),vcov(lm_interest_direct_aux),Terms = c(16,9:11), df = summary(lm_interest_direct_aux)$df[2]) , silent = TRUE)
   h01 = try(wald.test(b = as.vector(lm_interest_direct_aux$m$getPars()),vcov(lm_interest_direct_aux),Terms = c(15,6:8), df = summary(lm_interest_direct_aux)$df[2]), silent = TRUE) 
   if(class(h03) %in% 'try-error' | class(h02) %in% 'try-error' |class(h01) %in% 'try-error'){next} else{
   #H03 results:   
   F_wald_h03 = h03$result$Ftest[1] 
   pval_wald_h03_F = h03$result$Ftest[4]
   Chi_wald_h03 = h03$result$chi2[1] 
   pval_wald_h03_Chi = h03$result$chi2[3]
   #H02 results: 
   F_wald_h02 = h02$result$Ftest[1] 
   pval_wald_h02_F = h02$result$Ftest[4]
   Chi_wald_h02 = h02$result$chi2[1] 
   pval_wald_h02_Chi = h02$result$chi2[3]
   #H01 results:   
   F_wald_h01 = h01$result$Ftest[1] 
   pval_wald_h01_F = h01$result$Ftest[4]
   Chi_wald_h01 = h01$result$chi2[1] 
   pval_wald_h01_Chi = h01$result$chi2[3]
   if(pval_wald_h02_F <= pval_wald_h03_F && pval_wald_h02_F <= pval_wald_h01_F){str_model_F = "LSTR2-ESTR"} else{str_model_F = "LSTR1"}
   if(pval_wald_h02_Chi <= pval_wald_h03_Chi && pval_wald_h02_Chi <= pval_wald_h01_Chi){str_model_Chi = "LSTR2-ESTR"} else{str_model_Chi = "LSTR1"}
   final_aux = rbind(F_test_statistic_lm, pval_F_lm, Chi_test_statistic_lm, pval_chi_lm,test_lm,str_model_F,str_model_Chi)
   colnames(final_aux) = names_interest[i]
   rownames(final_aux) = c("F_test","F_test_pval","Chi_test", "Chi_test_pval","LM result", "STR model F", "STR model Chi")
   LM_interest_state = cbind(LM_interest_state, final_aux)}
   }
}
write.csv(LM_interest_state, 'LM_interest_state.csv')


# --------------------- Estimation of STR models -------------------------------      
#Initialize and register parallel backend 
stopCluster(cl)
n_cores = detectCores()-1
registerDoParallel(cores = n_cores)
cl = makeCluster(n_cores)
registerDoSNOW(cl)

# ------- Inflation Equation ------------- #
# Selection of variables where LM test detected non-linearity
names_inflation_direct = colnames(LM_inflation_direct)
names_inflation_state = colnames(LM_inflation_state)
selected_variables_direct = names_inflation_direct[which(LM_inflation_direct[c("LM result"),] == "Nonlinear")]
selected_variables_state = names_inflation_state[which(LM_inflation_state[c("LM result"),] == "Nonlinear")]
length(which(LM_inflation_state[c("LM result"),] == "Nonlinear" & LM_inflation_state[c("STR model F"),] == "LSTR2-ESTR"))
selected_variables = c(selected_variables_direct,selected_variables_state)

data_str_inf_state = data_state_inflation[,selected_variables] 
data_str_inf_state = xts(data_str_inf_state, index_dates_data_state_final)
data_str_inf_model = as.data.frame(cbind(data_state_final$inflation, data_state_final$resid_inf, data_state_final$inflation.l1, data_state_final$output_gap, data_state_final$output_gap.l1, data_state_final$real_interest_gap.l2,
                                           data_state_final$rer_gap, data_state_final$inflation_ext.l2))
data_str_inf_model = xts(data_str_inf_model, index_dates_data_state_final)
data_str_inf = cbind.xts(data_str_inf_model, data_str_inf_state)
 
obs = nrow(data_str_inf)
data_str_inf = as.data.frame(data_str_inf)
names_str_inf = colnames(data_str_inf)
#Steps (for each nonlinear or state variable)
#1. Calculate Grid for Gamma and C parameters
#2. USe the Grid to compute the transition function (logistic/exponential)
#3. Add transition function to the model and estimate model by OLS - NLS
#4. Obtain SSR for each iteration, and choose the parameters that minimize it
#5. Use parameters in 4 as starting values to NLS estimate of STR


LSTR1_INF = list()
ptm <- proc.time()
#For Loop for estimates of LSTR1 function 
for(i in eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state))){
   #Declare state variable
  state_variable = data_str_inf[,i]  
  #Grid for Gamma and C parameters
  gamma_grid = seq(0.1*sd(state_variable[2:obs]), 32*sd(state_variable[2:obs]), length.out = 32) #Reescaled by transition variable 
  c_grid = seq(0.1, 0.95, 0.02) #From 10th percentile to 90th percentile
  combinations = as.data.frame(crossing(gamma = 1:length(gamma_grid), c = 1:length(c_grid)))
  dim(combinations)
  
  #LSTR 1 estimation
  ssr_linear_str = NULL
  parameters_linear_str = NULL
  get_init_param_i = function(x){
     #Combinations of parameters
     gamma_cond <- gamma_grid[combinations[x,1]]
     c_cond <-quantile(state_variable[2:obs], probs = c_grid[combinations[x,2]]) #State variable transition point
     # Compute transition function for the state variable
     f_lstr1 = 1/(1+exp(-(gamma_cond/sd(state_variable[2:obs]))*(state_variable-c_cond))) #This becomes a variable 
     # Estimate model by OLS (without rho term) to get initial parameters for NLS estimation
     ols_aux = lm(inflation ~ 1 + inflation.l1 + output_gap + output_gap.l1 + real_interest_gap.l2 + rer_gap + inflation_ext.l2 
                  +f_lstr1+eval(inflation.l1*f_lstr1) + eval(output_gap*f_lstr1) + eval(output_gap.l1*f_lstr1) + eval(real_interest_gap.l2*f_lstr1) + eval(rer_gap*f_lstr1) + eval(inflation_ext.l2*f_lstr1),
                  data = cbind(data_str_inf,f_lstr1))
     par_ols = coef(ols_aux)
     par_ols[is.na(par_ols)] <- 0
     # Add variable to the estimation of linear model (by NLS):
     #From x13 to x18 are the STR interactions of model variables with state variables
     #From x19 to x24 are lagged interactions with rho term
     
     linear_str = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4
                            +b5*x5+b6*x6 
                            -rho*b1*x7-rho*b2*x8-rho*b3*x9-rho*b4*x10
                            -rho*b5*x11-rho*b6*x12
                            +b7*x13 + b8*x14 + b9*x15 + b10*x16 + b11*x17 + b12*x18 + b13*x19
                            - rho*b7*x20 - rho*b8*x21 - rho*b9*x22 - rho*b10*x23 - rho*b11*x24 - rho*b12*x25-rho*b12*x26
                            , 
                            start = list(rho = 0.01, b0 = 0.01, b1 = 0.67, b2 = -0.11, b3= 0.38, b4= -0.09, b5= 0.084, b6 = 0.02, 
                                         b7 = par_ols[8], b8 = par_ols[9], b9 = par_ols[10], b10 = par_ols[11], b11 = par_ols[12], b12 = par_ols[13], b13 = par_ols[14]),
                            data = list(y = data_str_inf$inflation[2:obs], y_l = data_str_inf$inflation[1:obs-1], 
                                        x1 = data_str_inf$inflation.l1[2:obs], x2 = data_str_inf$output_gap[2:obs], x3 = data_str_inf$output_gap.l1[2:obs], 
                                        x4 = data_str_inf$real_interest_gap.l2[2:obs], x5 = data_str_inf$rer_gap[2:obs], x6 = data_str_inf$inflation_ext.l2[2:obs],
                                        x7 = data_str_inf$inflation.l1[1:obs-1], x8 = data_str_inf$output_gap[1:obs-1], x9 = data_str_inf$output_gap.l1[1:obs-1], 
                                        x10 = data_str_inf$real_interest_gap.l2[1:obs-1], x11 = data_str_inf$rer_gap[1:obs-1], x12 = data_str_inf$inflation_ext.l2[1:obs-1],
                                        x13 = f_lstr1[2:obs], x14 = data_str_inf$inflation.l1[2:obs]*f_lstr1[2:obs], x15 = data_str_inf$output_gap[2:obs]*f_lstr1[2:obs], x16 = data_str_inf$output_gap.l1[2:obs]*f_lstr1[2:obs], 
                                        x17 = data_str_inf$real_interest_gap.l2[2:obs]*f_lstr1[2:obs], x18 = data_str_inf$rer_gap[2:obs]*f_lstr1[2:obs], x19 = data_str_inf$inflation_ext.l2[2:obs]*f_lstr1[2:obs],
                                        x20 = f_lstr1[1:obs-1], x21 = data_str_inf$inflation.l1[1:obs-1]*f_lstr1[1:obs-1], x22 = data_str_inf$output_gap[1:obs-1]*f_lstr1[1:obs-1], x23 = data_str_inf$output_gap.l1[1:obs-1]*f_lstr1[1:obs-1], 
                                        x24 = data_str_inf$real_interest_gap.l2[1:obs-1]*f_lstr1[1:obs-1], x25 = data_str_inf$rer_gap[1:obs-1]*f_lstr1[1:obs-1], x26 = data_str_inf$inflation_ext.l2[1:obs-1]*f_lstr1[1:obs-1]
                                        
                            ), 
                            trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200), nls.control(minFactor = 1/4096, warnOnly = TRUE)), silent = TRUE)
     if(try(summary(linear_str), silent = TRUE)[2] %in% "try-error"||class(linear_str) %in% "try-error"){return(c(rep(NA, length(par_ols)+2)))} else {
        # Get the Residual Sum of Squares of the Regression
        ssr_aux = sum((predict(linear_str) - data_str_inf$inflation[2:obs]) ^ 2)  ## residual sum of squares
        #ssr_linear_str = rbind(ssr_linear_str, ssr_aux)
        # Get the model parameters
        parameters_aux = linear_str$m$getPars()
        parameters_ssr_linear_str = rbind(c(), c(parameters_aux, ssr_aux))
        return(parameters_ssr_linear_str)
     }
  }
  system.time({param_ssr_matrix = foreach(j = 1:nrow(combinations), .combine = rbind, .packages = c('parallel', 'minpack.lm')) %dopar% {get_init_param_i(j)}})
  param_ssr_matrix = na.omit(param_ssr_matrix)
  index_selection = which(param_ssr_matrix[,ncol(param_ssr_matrix)] == min(param_ssr_matrix[,ncol(param_ssr_matrix)])) 
  par_str_ini = c(param_ssr_matrix[index_selection,-ncol(param_ssr_matrix)],gamma_grid[combinations[index_selection,1]],quantile(state_variable,probs = c_grid[combinations[index_selection,2]]))
 
  # Estimate by NLS LM
    nls_str_LM = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4+b5*x5+b6*x6 
                         -rho*b1*x7-rho*b2*x8-rho*b3*x9-rho*b4*x10-rho*b5*x11-rho*b6*x12
                         +(1/(1+exp(-(gamma/sd(state_variable[2:obs]))*(s_var-c))))*(b7 +b8*x1 + b9*x2 + b10*x3 + b11*x4 + b12*x5 + b13*x6)
                         -rho*(1/(1+exp(-(gamma/sd(state_variable[2:obs]))*(s_var_l-c))))*(b7 +b8*x7 +b9*x8 +b10*x9 + b11*x10 + b12*x11 + b13*x12)
                         , 
                         data = list(y = data_str_inf$inflation[2:obs], y_l = data_str_inf$inflation[1:obs-1], 
                                     x1 = data_str_inf$inflation.l1[2:obs], x2 = data_str_inf$output_gap[2:obs], x3 = data_str_inf$output_gap.l1[2:obs], 
                                     x4 = data_str_inf$real_interest_gap.l2[2:obs], x5 = data_str_inf$rer_gap[2:obs], x6 = data_str_inf$inflation_ext.l2[2:obs],
                                     x7 = data_str_inf$inflation.l1[1:obs-1], x8 = data_str_inf$output_gap[1:obs-1], x9 = data_str_inf$output_gap.l1[1:obs-1], 
                                     x10 = data_str_inf$real_interest_gap.l2[1:obs-1], x11 = data_str_inf$rer_gap[1:obs-1], x12 = data_str_inf$inflation_ext.l2[1:obs-1],
                                     s_var = state_variable[2:obs], s_var_l = state_variable[1:obs-1]) 
                                     ,
                         start = list(rho = par_str_ini[1], b0 = par_str_ini[2], b1 = par_str_ini[3], b2 = par_str_ini[4], b3= par_str_ini[5], b4= par_str_ini[6], b5= par_str_ini[7], b6 = par_str_ini[8], 
                                      b7 = par_str_ini[9], b8 = par_str_ini[10], b9 = par_str_ini[11], b10 = par_str_ini[12], b11 = par_str_ini[13], b12 = par_str_ini[14], b13 = par_str_ini[15], gamma = par_str_ini[16], c = par_str_ini[17]),
                         
                         lower = c(0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,gamma_grid[1],quantile(state_variable,probs = c_grid[1])),
                         upper = c(Inf, Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,quantile(state_variable,probs = c_grid[length(c_grid)])),
                         control = nls.lm.control(maxiter = 200),nls.control(minFactor = 1/4096, warnOnly = TRUE), jac = NULL), silent = TRUE)
  
  if(try(summary(nls_str_LM),silent = TRUE)[2] %in% "try-error" || class(nls_str_LM) %in% "try-error"){next} else {
  convergence_message=if(nls_str_LM$convInfo$isConv == TRUE){convergence_message = "converged"} else {"no convergence"}
  #Thing we need to get from the model:
  # Name of the transition variable 
  name_state_variable = names_str_inf[i]
  # Parameters (linear and nonlinear part) 
  parameters_aux = nls_str_LM$m$getPars()
  #Std error of parameter estimates
  std_parameters_aux = summary(nls_str_LM)$coefficients[,2]
  #Transition Function (dated and ordered)
  f_lstr1_aux = 1/(1+exp(-(parameters_aux[length(parameters_aux)-1]/sd(state_variable[2:obs]))*(state_variable[2:obs]-parameters_aux[length(parameters_aux)])))
   #plot(f_lstr1_aux, type = 'l', col = 'blue')
   #plot(state_variable[2:obs],f_lstr1_aux, col = 'blue')  
  f_lstr1_aux_ordered = sort(f_lstr1_aux, decreasing = FALSE)
   #plot(f_lstr1_aux_ordered, col = 'red')
  #Predicted values: 
  predicted_str_aux = predict(nls_str_LM)
  #BIC (to compare against linear models) 
  BIC_str_aux = BIC(nls_str_LM)
  BIC_criteria = if(BIC_str_aux >= BIC_inf_nls){BIC_criteria = "Linear"} else {BIC_criteria = "No Linear"}
  #AIC (to compare against alternative STR especifications)
  AIC_str_aux = AIC(nls_str_LM)
  AIC_Criteria = if(AIC_str_aux >= AIC_inf_nls){AIC_Criteria = "Linear"} else {AIC_Criteria = "No Linear"}
  # Residuals of the model (with 0 at the beginning): 
  residuals_str_aux = resid(nls_str_LM)
  # SSR (necessary for evaluation stage)
  ssr_str_aux = sum((predict(nls_str_LM) - data_str_inf$inflation[2:obs])^2)
  SSR_criteria = if(ssr_str_aux >= ssr_1_inf_nls){SSR_criteria = "Linear"} else {SSR_criteria = "No Linear"}
  #Gradient for autocorrelation test 
  gradient_aux = nls_str_LM$m$gradient()
  #Test of restriction of nonlinear parameters = 0
  F_test_nlparam = ((ssr_1_inf_nls-ssr_str_aux)/(summary(inflation_eq_nls)$df[2] -summary(nls_str_LM)$df[2]))/(ssr_str_aux/(summary(nls_str_LM)$df[2]))
  pval_f_nlparam = pf(F_test_nlparam,summary(inflation_eq_nls)$df[2]-summary(nls_str_LM)$df[2],summary(nls_str_LM)$df[2],lower.tail = FALSE)
  #Test of autocorrelation of residuals 
    #1. Orthogonalize errors: 
    ortho_error_autocorr = lm(residuals_str_aux ~ 1 +gradient_aux)
    reisduals_ortho_str_aux = resid(ortho_error_autocorr)
    #2. Matrix of lagged residuales (vt):
    vt_l_residuals = cbind(reisduals_ortho_str_aux,dplyr::lag(reisduals_ortho_str_aux,1),dplyr::lag(reisduals_ortho_str_aux,2),dplyr::lag(reisduals_ortho_str_aux,3),
                            dplyr::lag(reisduals_ortho_str_aux,4),dplyr::lag(reisduals_ortho_str_aux,5),dplyr::lag(reisduals_ortho_str_aux,6),dplyr::lag(reisduals_ortho_str_aux,7),
                            dplyr::lag(reisduals_ortho_str_aux,9),dplyr::lag(reisduals_ortho_str_aux,10),dplyr::lag(reisduals_ortho_str_aux,11),dplyr::lag(reisduals_ortho_str_aux,12))
    vt_l_residuals[is.na(vt_l_residuals)] <- 0
    vt_l_residuals = vt_l_residuals[,c(-1)]
    #3. Regression result and LM tests (F and Chi) 
    auto_corr_test = lm(residuals_str_aux ~ gradient_aux + vt_l_residuals)
    ssr_alt = sum((predict(auto_corr_test)-residuals_str_aux)^2) 
    F_corr_str = ((ssr_str_aux - ssr_alt)/(summary(nls_str_LM)$df[2]-auto_corr_test$df.residual))/(ssr_alt/(auto_corr_test$df.residual))
    pval_f_corr = pf(F_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual,auto_corr_test$df.residual,lower.tail = FALSE)
    chi_corr_str = length(residuals_str_aux)*((ssr_str_aux - ssr_alt)/ssr_str_aux) 
    pval_chi_corr = pchisq(chi_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual, lower.tail = FALSE) 
  #Result of Jarque-Bera and Shapiro Test 
  jb_test_aux = jarque.bera.test(residuals_str_aux)$p.value
  shapiro_aux = shapiro.test(residuals_str_aux)$p.value
  #Ramsey Test
  ols_str_ramsey = lm(residuals_str_aux ~ gradient_aux + predicted_str_aux + predicted_str_aux^2 + predicted_str_aux^3)
  ssr_ols_str_ramsey = sum((predict(ols_str_ramsey) - residuals_str_aux)^2)
    #F test statistic:
    f_ramsey_ols = ((ssr_str_aux - ssr_ols_str_ramsey)/(summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual))/(ssr_ols_str_ramsey/(ols_str_ramsey$df.residual))
    pval_f_str_ramsey = pf(f_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,ols_str_ramsey$df.residual, lower.tail = FALSE)
    #Chi test statistic 
    chi_ramsey_ols = length(residuals_str_aux)*(ssr_str_aux - ssr_ols_str_ramsey)/ssr_str_aux
    pval_chi_str_ramsey = pchisq(chi_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,lower.tail = FALSE)
  #No additive Linearity Test - LM: 
   no_ad_nolin_results_direct = NULL
   for(k in eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct))){
   state_2 =  data_str_inf[2:nrow(data_str_inf),k]
   no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux + eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
   ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
   f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
   pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
   no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
   no_ad_nolin_results_direct = rbind(no_ad_nolin_results_direct,no_ad_nolin_results_aux)
   }
   
   no_ad_nolin_results_state = NULL
   for(k in eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state))){
   state_2 =  data_str_inf[2:nrow(data_str_inf),k]
   no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux +state_2+eval(state_2^2)+eval(state_2^3)+ eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
   ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
   f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
   pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
   no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
   no_ad_nolin_results_state = rbind(no_ad_nolin_results_state,no_ad_nolin_results_aux)
   }
   no_ad_nolin_results_direct = na.omit(no_ad_nolin_results_direct)
   no_ad_nolin_results_state = na.omit(no_ad_nolin_results_state)
  if(any(no_ad_nolin_results_direct[,2] <= 0.05) == TRUE){no_ad_d_nl_5_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_5_pc="No Additive Nonlinearity"}
  if(any(no_ad_nolin_results_direct[,2] <= 0.01) == TRUE){no_ad_d_nl_1_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_1_pc="No Additive Nonlinearity"}
  if(any(no_ad_nolin_results_state[,2] <= 0.05) == TRUE){no_ad_s_nl_5_pc = length(which(no_ad_nolin_results_state[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_5_pc="No Additive Nonlinearity"}
  if(any(no_ad_nolin_results_state[,2] <= 0.01) == TRUE){no_ad_s_nl_1_pc = length(which(no_ad_nolin_results_state[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_1_pc="No Additive Nonlinearity"}
  #Parameter Constancy Test - LM: cannot be performed do to small sample size problem (low DF)
  #Arch test: 
  noarch_test =  ArchTest(residuals_str_aux, lags = 12)
  arch_chi_statistic = noarch_test$statistic
  arch_chi_pval = noarch_test$p.value
  
 #Final Output
 model_info = rbind(convergence_message,BIC_str_aux,BIC_criteria, AIC_str_aux, AIC_Criteria, ssr_str_aux,SSR_criteria,F_corr_str,pval_f_corr,chi_corr_str,pval_chi_corr,
                   jb_test_aux,shapiro_aux,f_ramsey_ols,pval_f_str_ramsey,chi_ramsey_ols,pval_chi_str_ramsey,no_ad_d_nl_5_pc,no_ad_d_nl_1_pc,no_ad_s_nl_5_pc,no_ad_s_nl_1_pc,arch_chi_statistic,arch_chi_pval, 
                   F_test_nlparam, pval_f_nlparam) 
 colnames(model_info) = name_state_variable
 rownames(model_info) = c("Convergence", "BIC", "BIC Criteria", "AIC", "AIC Criteria", "SSR", "SSR Criteria", "F test corr", "Pval F test", "Chi2 test corr", "Pval Chi2 test",
                          "JB test", "Shapiro test", "F test Ramsey", "Pval F Ramsey", "Chi2 test Ramsey", "Pval Chi2 Ramsey", "No Ad NL D 5%","No Ad D NL 1%","No Ad NL S 5%","No Ad S NL 1%", "ARCH test stat", "ARCH test pval", 
                          "F test nl", "Pval F nl")
 Transition_function = cbind(f_lstr1_aux, f_lstr1_aux_ordered, state_variable[2:obs])
 colnames(Transition_function) = c("LSTR1", "LSTR1 sorted", "State_Variable")
 final_output = list(Model_Info = model_info,Model = nls_str_LM, Transition_function = Transition_function)
 LSTR1_INF[[eval(i-ncol(data_str_inf_model))]] <- final_output
 }
}
names(LSTR1_INF) <- names_str_inf[eval(ncol(data_str_inf_model)+1):length(names_str_inf)]
proc.time() - ptm 

#Results: of model info 
LSTR1_INF_model_results = NULL
for(j in 1:length(LSTR1_INF)){aux_info = LSTR1_INF[[j]]$Model_Info
LSTR1_INF_model_results = cbind(LSTR1_INF_model_results,aux_info)
}
#Selected models According to BIC/AIC criteria,  
LSTR1_INF_selected_models = LSTR1_INF_model_results[,
which(LSTR1_INF_model_results['Convergence',] == "converged" & 
         (LSTR1_INF_model_results['BIC Criteria',] == "No Linear" |  LSTR1_INF_model_results['AIC Criteria',] == "No Linear")
      & (LSTR1_INF_model_results['JB test',] >= 0.05 | LSTR1_INF_model_results['Shapiro test',] >= 0.05)
      &  LSTR1_INF_model_results['Pval F nl',] <= 0.05 &  LSTR1_INF_model_results['Pval F test',] >= 0.05 &  LSTR1_INF_model_results['Pval Chi2 test',] >= 0.05
      &  LSTR1_INF_model_results['Pval F Ramsey',] >= 0.05 & LSTR1_INF_model_results['Pval Chi2 Ramsey',] >= 0.05
      & LSTR1_INF_model_results['ARCH test pval',] >= 0.05 & (LSTR1_INF_model_results['No Ad D NL 1%',] == "No Additive Nonlinearity" | LSTR1_INF_model_results['No Ad D NL 1%',] < 0.25) 
      & (LSTR1_INF_model_results['No Ad S NL 1%',] == "No Additive Nonlinearity" | LSTR1_INF_model_results['No Ad S NL 1%',] < 0.25))]
write.csv(LSTR1_INF_selected_models, "LSTR1_INF_selected.csv")
# Selection of variables where LM test detected non-linearity and Wald test suggested LSTR2 - ESTR model estimation
names_inflation_direct = colnames(LM_inflation_direct)
names_inflation_state = colnames(LM_inflation_state)
selected_variables_direct = names_inflation_direct[which(LM_inflation_direct[c("LM result"),] == "Nonlinear" & LM_inflation_direct[c("STR model F"),] == "LSTR2-ESTR")]
selected_variables_state = names_inflation_state[which(LM_inflation_state[c("LM result"),] == "Nonlinear" & LM_inflation_state[c("STR model F"),] == "LSTR2-ESTR")]
length(which(LM_inflation_state[c("LM result"),] == "Nonlinear" & LM_inflation_state[c("STR model F"),] == "LSTR2-ESTR"))
selected_variables = c(selected_variables_direct,selected_variables_state)

data_str_inf_state = data_state_inflation[,selected_variables] 
data_str_inf_state = xts(data_str_inf_state, index_dates_data_state_final)
data_str_inf_model = as.data.frame(cbind(data_state_final$inflation, data_state_final$resid_inf, data_state_final$inflation.l1, data_state_final$output_gap, data_state_final$output_gap.l1, data_state_final$real_interest_gap.l2,
                                         data_state_final$rer_gap, data_state_final$inflation_ext.l2))
data_str_inf_model = xts(data_str_inf_model, index_dates_data_state_final)
data_str_inf = cbind.xts(data_str_inf_model, data_str_inf_state)

obs = nrow(data_str_inf)
data_str_inf = as.data.frame(data_str_inf)
names_str_inf = colnames(data_str_inf)


ESTR_INF = list()
ptm <- proc.time()
#For Loop for estimates of LSTR1 function 
for(i in eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state))){
   #Declare state variable
   state_variable = data_str_inf[,i]  
   #Grid for Gamma and C parameters
      gamma_grid = seq(0.1*sd(state_variable[2:obs]), 32*sd(state_variable[2:obs]), length.out = 32) #Reescaled by transition variable 
      c_grid = seq(0.1, 0.95, 0.02) #From 10th percentile to 90th percentile
   combinations = as.data.frame(crossing(gamma = 1:length(gamma_grid), c = 1:length(c_grid)))
   dim(combinations)
   #ESTR  estimation
   ssr_linear_str = NULL
   parameters_linear_str = NULL
   get_init_param_i = function(x){
      #Combinations of parameters
      gamma_cond <- gamma_grid[combinations[x,1]]
      c_cond <-quantile(state_variable, probs = c_grid[combinations[x,2]]) #State variable transition point
      # Compute transition function for the state variable
      f_estr1 = 1-exp(-(gamma_cond/sd(state_variable[2:obs]))*((state_variable-c_cond)^2)) #This becomes a variable 
      # Estimate model by OLS (without rho term) to get initial parameters for NLS estimation
      ols_aux = lm(inflation ~ 1 + inflation.l1 + output_gap + output_gap.l1 + real_interest_gap.l2 + rer_gap + inflation_ext.l2 
                   +f_estr1+eval(inflation.l1*f_estr1) + eval(output_gap*f_estr1) + eval(output_gap.l1*f_estr1) + eval(real_interest_gap.l2*f_estr1) + eval(rer_gap*f_estr1) + eval(inflation_ext.l2*f_estr1),
                   data = data_str_inf)
      par_ols = coef(ols_aux)
      par_ols[is.na(par_ols)] <- 0
      # Add variable to the estimation of linear model (by NLS):
      #From x13 to x18 are the STR interactions of model variables with state variables
      #From x19 to x24 are lagged interctions with rho term
      #b7 = par_ols[8], b8 = par_ols[9], b9 = par_ols[10], b10 = par_ols[11], b11 = par_ols[12], b12 = par_ols[13], b13 = par_ols[14]),
      
      linear_str = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4
                             +b5*x5+b6*x6 
                             -rho*b1*x7-rho*b2*x8-rho*b3*x9-rho*b4*x10
                             -rho*b5*x11-rho*b6*x12
                             +b7*x13 + b8*x14 + b9*x15 + b10*x16 + b11*x17 + b12*x18 + b13*x19
                             - rho*b7*x20 - rho*b8*x21 - rho*b9*x22 - rho*b10*x23 - rho*b11*x24 - rho*b12*x25 - rho*b13*x26
                             , 
                             start = list(rho = 0.01, b0 = 0.01, b1 = 0.67, b2 = -0.11, b3= 0.38, b4= -0.09, b5= 0.084, b6 = 0.02, 
                                          b7 = par_ols[8], b8 = par_ols[9], b9 = par_ols[10], b10 = par_ols[11], b11 = par_ols[12], b12 = par_ols[13], b13 = par_ols[14]),
                             data = list(y = data_str_inf$inflation[2:obs], y_l = data_str_inf$inflation[1:obs-1], 
                                         x1 = data_str_inf$inflation.l1[2:obs], x2 = data_str_inf$output_gap[2:obs], x3 = data_str_inf$output_gap.l1[2:obs], 
                                         x4 = data_str_inf$real_interest_gap.l2[2:obs], x5 = data_str_inf$rer_gap[2:obs], x6 = data_str_inf$inflation_ext.l2[2:obs],
                                         x7 = data_str_inf$inflation.l1[1:obs-1], x8 = data_str_inf$output_gap[1:obs-1], x9 = data_str_inf$output_gap.l1[1:obs-1], 
                                         x10 = data_str_inf$real_interest_gap.l2[1:obs-1], x11 = data_str_inf$rer_gap[1:obs-1], x12 = data_str_inf$inflation_ext.l2[1:obs-1],
                                         x13 = f_estr1[2:obs],x14 = data_str_inf$inflation.l1[2:obs]*f_estr1[2:obs], x15 = data_str_inf$output_gap[2:obs]*f_estr1[2:obs], x16 = data_str_inf$output_gap.l1[2:obs]*f_estr1[2:obs], 
                                         x17 = data_str_inf$real_interest_gap.l2[2:obs]*f_estr1[2:obs], x18 = data_str_inf$rer_gap[2:obs]*f_estr1[2:obs], x19 = data_str_inf$inflation_ext.l2[2:obs]*f_estr1[2:obs],
                                         x20 = f_estr1[1:obs-1], x21 = data_str_inf$inflation.l1[1:obs-1]*f_estr1[1:obs-1], x22 = data_str_inf$output_gap[1:obs-1]*f_estr1[1:obs-1], x23 = data_str_inf$output_gap.l1[1:obs-1]*f_estr1[1:obs-1], 
                                         x24 = data_str_inf$real_interest_gap.l2[1:obs-1]*f_estr1[1:obs-1], x25 = data_str_inf$rer_gap[1:obs-1]*f_estr1[1:obs-1], x26 = data_str_inf$inflation_ext.l2[1:obs-1]*f_estr1[1:obs-1]
                                         
                             ), 
                             trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200), nls.control(minFactor = 1/4096, warnOnly = TRUE)), silent = TRUE)
      if(try(summary(linear_str), silent = TRUE)[2] %in% "try-error" ||class(linear_str) %in% "try-error"){return(c(rep(NA, length(par_ols)+2)))} else {
         # Get the Residual Sum of Squares of the Regression
         ssr_aux = sum((predict(linear_str) - data_str_inf$inflation[2:obs]) ^ 2)  ## residual sum of squares
         #ssr_linear_str = rbind(ssr_linear_str, ssr_aux)
         # Get the model parameters
         parameters_aux = linear_str$m$getPars()
         parameters_ssr_linear_str = rbind(c(), c(parameters_aux, ssr_aux))
         return(parameters_ssr_linear_str)
      }  
      
   }
   system.time({param_ssr_matrix = foreach(i = 1:nrow(combinations), .combine = rbind,.packages = c('doParallel', 'minpack.lm')) %dopar% {get_init_param_i(i)}})
   param_ssr_matrix = na.omit(param_ssr_matrix)
   index_selection = which(param_ssr_matrix[,ncol(param_ssr_matrix)] == min(param_ssr_matrix[,ncol(param_ssr_matrix)])) 
   par_str_ini = c(param_ssr_matrix[index_selection,-ncol(param_ssr_matrix)],gamma_grid[combinations[index_selection,1]],quantile(state_variable,probs = c_grid[combinations[index_selection,2]]))
   # Estimate by NLS LM
   nls_str_LM = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4
                              +b5*x5+b6*x6 
                              -rho*b1*x7-rho*b2*x8-rho*b3*x9-rho*b4*x10
                              -rho*b5*x11-rho*b6*x12
                              +(1-exp(-(gamma/sd(state_variable[2:obs]))*((s_var-c)^2)))*(b7 + b8*x1 + b9*x2 + b10*x3 + b11*x4 + b12*x5 + b13*x6)
                              -rho*(1-exp(-(gamma/sd(state_variable[2:obs]))*((s_var_l-c)^2)))*(b7 +b8*x7 + b9*x8 + b10*x9 + b11*x10 + b12*x11 + b13*x12)
                              , 
                              data = list(y = data_str_inf$inflation[2:obs], y_l = data_str_inf$inflation[1:obs-1], 
                                          x1 = data_str_inf$inflation.l1[2:obs], x2 = data_str_inf$output_gap[2:obs], x3 = data_str_inf$output_gap.l1[2:obs], 
                                          x4 = data_str_inf$real_interest_gap.l2[2:obs], x5 = data_str_inf$rer_gap[2:obs], x6 = data_str_inf$inflation_ext.l2[2:obs],
                                          x7 = data_str_inf$inflation.l1[1:obs-1], x8 = data_str_inf$output_gap[1:obs-1], x9 = data_str_inf$output_gap.l1[1:obs-1], 
                                          x10 = data_str_inf$real_interest_gap.l2[1:obs-1], x11 = data_str_inf$rer_gap[1:obs-1], x12 = data_str_inf$inflation_ext.l2[1:obs-1],
                                          s_var = state_variable[2:obs], s_var_l = state_variable[1:obs-1]) 
                              ,
                              start = list(rho = par_str_ini[1], b0 = par_str_ini[2], b1 = par_str_ini[3], b2 = par_str_ini[4], b3= par_str_ini[5], b4= par_str_ini[6], b5= par_str_ini[7], b6 = par_str_ini[8], 
                                           b7 = par_str_ini[9], b8 = par_str_ini[10], b9 = par_str_ini[11], b10 = par_str_ini[12], b11 = par_str_ini[13], b12 = par_str_ini[14],b13 = par_str_ini[15], gamma = par_str_ini[16], c = par_str_ini[17]),
                              
                              lower = c(0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf, - Inf, -Inf,-Inf,gamma_grid[1],quantile(state_variable,probs = c_grid[1])),
                              upper = c(Inf, Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf, quantile(state_variable,probs = c_grid[length(c_grid)])),
                              control = nls.lm.control(maxiter = 200), nls.control(minFactor = 1/4096, warnOnly = TRUE), jac = NULL), silent = TRUE)
    
   if(try(summary(nls_str_LM),silent = TRUE)[2] %in% "try-error" || class(nls_str_LM) %in% "try-error"){next} else {
      convergence_message=if(nls_str_LM$convInfo$isConv == TRUE){convergence_message = "converged"} else {"no convergence"}
      #Thing we need to get from the model:
      # Name of the transition variable 
      name_state_variable = names_str_inf[i]
      # Parameters (linear and nonlinear part) 
      parameters_aux = nls_str_LM$m$getPars()
      #Std error of parameter estimates
      std_parameters_aux = summary(nls_str_LM)$coefficients[,2]
      #Transition Function (dated and ordered)
      f_estr1_aux = (1-exp(-(parameters_aux[length(parameters_aux)-1]/sd(state_variable[2:obs]))*(state_variable[2:obs]-parameters_aux[length(parameters_aux)])^2))
         #plot(f_estr1_aux, type = 'l', col = 'blue')
         #plot(state_variable[2:obs],f_estr1_aux, col = 'blue')
      f_estr1_aux_ordered = sort(f_estr1_aux, decreasing = FALSE)
         #plot(f_estr1_aux_ordered, col = 'red')
      #Predicted values: 
      predicted_str_aux = predict(nls_str_LM)
      #BIC (to compare against linear models) 
      BIC_str_aux = BIC(nls_str_LM)
      BIC_criteria = if(BIC_str_aux >= BIC_inf_nls){BIC_criteria = "Linear"} else {BIC_criteria = "No Linear"}
      #AIC (to compare against alternative STR especifications)
      AIC_str_aux = AIC(nls_str_LM)
      AIC_Criteria = if(AIC_str_aux >= AIC_inf_nls){AIC_Criteria = "Linear"} else {AIC_Criteria = "No Linear"}
      # Residuals of the model (with 0 at the beginning): 
      residuals_str_aux = resid(nls_str_LM)
      # SSR (necessary for evaluation stage)
      ssr_str_aux = sum((predict(nls_str_LM) - data_str_inf$inflation[2:obs])^2)
      SSR_criteria = if(ssr_str_aux >= ssr_1_inf_nls){SSR_criteria = "Linear"} else {SSR_criteria = "No Linear"}
      #Gradient for autocorrelation test 
      gradient_aux = nls_str_LM$m$gradient()
      #Test of restriction of nonlinear parameters = 0
      F_test_nlparam = ((ssr_1_inf_nls-ssr_str_aux)/(summary(inflation_eq_nls)$df[2] -summary(nls_str_LM)$df[2]))/(ssr_str_aux/(summary(nls_str_LM)$df[2]))
      pval_f_nlparam = pf(F_test_nlparam,summary(inflation_eq_nls)$df[2]-summary(nls_str_LM)$df[2],summary(nls_str_LM)$df[2],lower.tail = FALSE)
      #Test of autocorrelation of residuals 
      #1. Orthogonalize errors: 
      ortho_error_autocorr = lm(residuals_str_aux ~ 1 +gradient_aux)
      reisduals_ortho_str_aux = resid(ortho_error_autocorr)
      #2. Matrix of lagged residuales (vt):
      vt_l_residuals = cbind(reisduals_ortho_str_aux,dplyr::lag(reisduals_ortho_str_aux,1),dplyr::lag(reisduals_ortho_str_aux,2),dplyr::lag(reisduals_ortho_str_aux,3),
                             dplyr::lag(reisduals_ortho_str_aux,4),dplyr::lag(reisduals_ortho_str_aux,5),dplyr::lag(reisduals_ortho_str_aux,6),dplyr::lag(reisduals_ortho_str_aux,7),
                             dplyr::lag(reisduals_ortho_str_aux,9),dplyr::lag(reisduals_ortho_str_aux,10),dplyr::lag(reisduals_ortho_str_aux,11),dplyr::lag(reisduals_ortho_str_aux,12))
      vt_l_residuals[is.na(vt_l_residuals)] <- 0
      vt_l_residuals = vt_l_residuals[,c(-1)]
      #3. Regression result and LM tests (F and Chi) 
      auto_corr_test = lm(residuals_str_aux ~ gradient_aux + vt_l_residuals)
      ssr_alt = sum((predict(auto_corr_test)-residuals_str_aux)^2) 
      F_corr_str = ((ssr_str_aux - ssr_alt)/(summary(nls_str_LM)$df[2]-auto_corr_test$df.residual))/(ssr_alt/(auto_corr_test$df.residual))
      pval_f_corr = pf(F_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual,auto_corr_test$df.residual,lower.tail = FALSE)
      chi_corr_str = length(residuals_str_aux)*((ssr_str_aux - ssr_alt)/ssr_str_aux) 
      pval_chi_corr = pchisq(chi_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual, lower.tail = FALSE) 
      #Result of Jarque-Bera and Shapiro Test 
      jb_test_aux = jarque.bera.test(residuals_str_aux)$p.value
      shapiro_aux = shapiro.test(residuals_str_aux)$p.value
      #Ramsey Test
      ols_str_ramsey = lm(residuals_str_aux ~ gradient_aux + predicted_str_aux + predicted_str_aux^2 + predicted_str_aux^3)
      ssr_ols_str_ramsey = sum((predict(ols_str_ramsey) - residuals_str_aux)^2)
      #F test statistic:
      f_ramsey_ols = ((ssr_str_aux - ssr_ols_str_ramsey)/(summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual))/(ssr_ols_str_ramsey/(ols_str_ramsey$df.residual))
      pval_f_str_ramsey = pf(f_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,ols_str_ramsey$df.residual, lower.tail = FALSE)
      #Chi test statistic 
      chi_ramsey_ols = length(residuals_str_aux)*(ssr_str_aux - ssr_ols_str_ramsey)/ssr_str_aux
      pval_chi_str_ramsey = pchisq(chi_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,lower.tail = FALSE)
      #No additive Linearity Test - LM: 
      no_ad_nolin_results_direct = NULL
      for(k in eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct))){
         state_2 =  data_str_inf[2:nrow(data_str_inf),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux + eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_direct = rbind(no_ad_nolin_results_direct,no_ad_nolin_results_aux)
      }
      
      no_ad_nolin_results_state = NULL
      for(k in eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state))){
         state_2 =  data_str_inf[2:nrow(data_str_inf),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux +state_2+eval(state_2^2)+eval(state_2^3)+ eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_state = rbind(no_ad_nolin_results_state,no_ad_nolin_results_aux)
      }
      no_ad_nolin_results_direct = na.omit(no_ad_nolin_results_direct)
      no_ad_nolin_results_state = na.omit(no_ad_nolin_results_state)
      if(any(no_ad_nolin_results_direct[,2] <= 0.05) == TRUE){no_ad_d_nl_5_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_direct[,2] <= 0.01) == TRUE){no_ad_d_nl_1_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_1_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.05) == TRUE){no_ad_s_nl_5_pc = length(which(no_ad_nolin_results_state[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.01) == TRUE){no_ad_s_nl_1_pc = length(which(no_ad_nolin_results_state[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_1_pc="No Additive Nonlinearity"}
      #Parameter Constancy Test - LM: cannot be performed do to small sample size problem (low DF)
      #Arch test: 
      noarch_test =  ArchTest(residuals_str_aux, lags = 12)
      arch_chi_statistic = noarch_test$statistic
      arch_chi_pval = noarch_test$p.value
      
      #Final Output
      model_info = rbind(convergence_message,BIC_str_aux,BIC_criteria, AIC_str_aux, AIC_Criteria, ssr_str_aux,SSR_criteria,F_corr_str,pval_f_corr,chi_corr_str,pval_chi_corr,
                         jb_test_aux,shapiro_aux,f_ramsey_ols,pval_f_str_ramsey,chi_ramsey_ols,pval_chi_str_ramsey,no_ad_d_nl_5_pc,no_ad_d_nl_1_pc,no_ad_s_nl_5_pc,no_ad_s_nl_1_pc,arch_chi_statistic,arch_chi_pval, 
                         F_test_nlparam, pval_f_nlparam)
       
      colnames(model_info) = name_state_variable
      rownames(model_info) = c("Convergence", "BIC", "BIC Criteria", "AIC", "AIC Criteria", "SSR", "SSR Criteria", "F test corr", "Pval F test", "Chi2 test corr", "Pval Chi2 test",
                               "JB test", "Shapiro test", "F test Ramsey", "Pval F Ramsey", "Chi2 test Ramsey", "Pval Chi2 Ramsey", "No Ad NL D 5%","No Ad D NL 1%","No Ad NL S 5%","No Ad S NL 1%", "ARCH test stat", "ARCH test pval", 
                               "F test nl", "Pval F nl")
      Transition_function = cbind(f_estr1_aux, f_estr1_aux_ordered)
      colnames(Transition_function) = c("LSTR1", "LSTR1 sorted")
      final_output = list(Model_Info = model_info,Model = nls_str_LM, Transition_function = Transition_function)
      ESTR_INF[[eval(i-ncol(data_str_inf_model))]] <- final_output
   }
}
names(ESTR_INF) <- names_str_inf[eval(ncol(data_str_inf_model)+1):length(names_str_inf)]
proc.time() - ptm 

#Results: of model info 
ESTR_INF_model_results = NULL
for(j in 1:length(ESTR_INF)){aux_info = ESTR_INF[[j]]$Model_Info
ESTR_INF_model_results = cbind(ESTR_INF_model_results,aux_info)
}
#Selected models According to BIC/AIC criteria,  
ESTR_INF_selected_models = ESTR_INF_model_results[,
                                                    which(ESTR_INF_model_results['Convergence',] == "converged" & 
                                                             (ESTR_INF_model_results['BIC Criteria',] == "No Linear" |  ESTR_INF_model_results['AIC Criteria',] == "No Linear")
                                                          & (ESTR_INF_model_results['JB test',] >= 0.05 | ESTR_INF_model_results['Shapiro test',] >= 0.05)
                                                          &  ESTR_INF_model_results['Pval F nl',] <= 0.05 &  ESTR_INF_model_results['Pval F test',] >= 0.05 &  ESTR_INF_model_results['Pval Chi2 test',] >= 0.05
                                                          &  ESTR_INF_model_results['Pval F Ramsey',] >= 0.05 & ESTR_INF_model_results['Pval Chi2 Ramsey',] >= 0.05
                                                          & ESTR_INF_model_results['ARCH test pval',] >= 0.05 & (ESTR_INF_model_results['No Ad D NL 1%',] == "No Additive Nonlinearity" | ESTR_INF_model_results['No Ad D NL 1%',] < 0.25) 
                                                          & (ESTR_INF_model_results['No Ad S NL 1%',] == "No Additive Nonlinearity" | ESTR_INF_model_results['No Ad S NL 1%',] < 0.25))]
 
write.csv(ESTR_INF_selected_models, "ESTR_INF_selected.csv")

LSTR2_INF = list()
ptm <- proc.time()
#For Loop for estimates of LSTR1 function 
for(i in eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state))){
   #Declare state variable
   state_variable = data_str_inf[,i]  
   #Grid for Gamma and C parameters
   gamma_grid = seq(0.1*sd(state_variable[2:obs])^2, 32*sd(state_variable[2:obs])^2, length.out = 32) #Reescaled by transition variable 
   c1_grid = seq(0.1,0.5,0.04) #c1 parameter, from 10 th percentile to 50th percentile 
   c2_grid = seq(0.55, 0.95, 0.04) #c2 parameter, from 50th percentile to 90th percentile 
   combinations = as.data.frame(crossing(gamma = 1:length(gamma_grid), c1 = 1:length(c1_grid), c2 = 1:length(c2_grid)))
   dim(combinations)
   
   #LSTR 2 estimation
   ssr_linear_str = NULL
   parameters_linear_str = NULL
   get_init_param_i <- function(x){
      #Combinations of parameters
      gamma_cond <- gamma_grid[combinations[x,1]]
      c1_cond <-quantile(state_variable, probs = c1_grid[combinations[x,2]]) #State variable transition point
      c2_cond <-quantile(state_variable, probs = c2_grid[combinations[x,3]]) #State variable transition point
      # Compute transition function for the state variable
      f_lstr2 = 1/(1+exp(-(gamma_cond/sd(state_variable[2:obs])^2)*(state_variable-c1_cond)*(state_variable-c2_cond))) #This becomes a variable 
      # Estimate model by OLS (without rho term) to get initial parameters for NLS estimation
      ols_aux = lm(inflation ~ 1 + inflation.l1 + output_gap + output_gap.l1 + real_interest_gap.l2 + rer_gap + inflation_ext.l2 
                   +f_lstr2+ eval(inflation.l1*f_lstr2) + eval(output_gap*f_lstr2) + eval(output_gap.l1*f_lstr2) + eval(real_interest_gap.l2*f_lstr2) + eval(rer_gap*f_lstr2) + eval(inflation_ext.l2*f_lstr2),
                   data = data_str_inf)
      par_ols = coef(ols_aux)
      par_ols[is.na(par_ols)]<- 0
      # Add variable to the estimation of linear model (by NLS):
      #From x13 to x18 are the STR interactions of model variables with state variables
      #From x19 to x24 are lagged interctions with rho term
      #b7 = par_ols[8], b8 = par_ols[9], b9 = par_ols[10], b10 = par_ols[11], b11 = par_ols[12], b12 = par_ols[13], b13 = par_ols[14]),
      
      linear_str = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4
                             +b5*x5+b6*x6 
                             -rho*b1*x7-rho*b2*x8-rho*b3*x9-rho*b4*x10
                             -rho*b5*x11-rho*b6*x12
                             +b7*x13 + b8*x14 + b9*x15 + b10*x16 + b11*x17 + b12*x18 + b13*x19
                             - rho*b7*x20 - rho*b8*x21 - rho*b9*x22 - rho*b10*x23 - rho*b11*x24 - rho*b12*x25 - rho*b13*x26
                             , 
                             start = list(rho = 0.01, b0 = 0.01, b1 = 0.67, b2 = -0.11, b3= 0.38, b4= -0.09, b5= 0.084, b6 = 0.02, 
                                          b7 = par_ols[8], b8 = par_ols[9], b9 = par_ols[10], b10 = par_ols[11], b11 = par_ols[12], b12 = par_ols[13], b13 = par_ols[14]),
                             data = list(y = data_str_inf$inflation[2:obs], y_l = data_str_inf$inflation[1:obs-1], 
                                         x1 = data_str_inf$inflation.l1[2:obs], x2 = data_str_inf$output_gap[2:obs], x3 = data_str_inf$output_gap.l1[2:obs], 
                                         x4 = data_str_inf$real_interest_gap.l2[2:obs], x5 = data_str_inf$rer_gap[2:obs], x6 = data_str_inf$inflation_ext.l2[2:obs],
                                         x7 = data_str_inf$inflation.l1[1:obs-1], x8 = data_str_inf$output_gap[1:obs-1], x9 = data_str_inf$output_gap.l1[1:obs-1], 
                                         x10 = data_str_inf$real_interest_gap.l2[1:obs-1], x11 = data_str_inf$rer_gap[1:obs-1], x12 = data_str_inf$inflation_ext.l2[1:obs-1],
                                         x13 = f_lstr2[2:obs], x14 = data_str_inf$inflation.l1[2:obs]*f_lstr2[2:obs], x15 = data_str_inf$output_gap[2:obs]*f_lstr2[2:obs], x16 = data_str_inf$output_gap.l1[2:obs]*f_lstr2[2:obs], 
                                         x17 = data_str_inf$real_interest_gap.l2[2:obs]*f_lstr2, x18 = data_str_inf$rer_gap[2:obs]*f_lstr2, x19 = data_str_inf$inflation_ext.l2[2:obs]*f_lstr2[2:obs],
                                         x20 = f_lstr2[1:obs-1], x21 = data_str_inf$inflation.l1[1:obs-1]*f_lstr2[1:obs-1], x22 = data_str_inf$output_gap[1:obs-1]*f_lstr2[1:obs-1], x23 = data_str_inf$output_gap.l1[1:obs-1]*f_lstr2[1:obs-1], 
                                         x24 = data_str_inf$real_interest_gap.l2[1:obs-1]*f_lstr2[1:obs-1], x25 = data_str_inf$rer_gap[1:obs-1]*f_lstr2[1:obs-1], x26 = data_str_inf$inflation_ext.l2[1:obs-1]*f_lstr2[1:obs-1]
                                         
                             ), 
                             control = nls.lm.control(maxiter = 200,), jac = NULL, trace =  TRUE, nls.control(minFactor = 1/4096, warnOnly = TRUE)), silent = TRUE)
      if(try(summary(linear_str), silent = TRUE)[2] %in% "try-error" ||class(linear_str) %in% "try-error"){return(rep(NA, 16))} else {
         # Get the Residual Sum of Squares of the Regression
         ssr_aux = sum((predict(linear_str) - data_str_inf$inflation[2:obs])^2)  ## residual sum of squares
         #ssr_linear_str = rbind(ssr_linear_str, ssr_aux)
         # Get the model parameters
         parameters_aux = linear_str$m$getPars()
         parameters_ssr_linear_str = rbind(c(), c(parameters_aux, ssr_aux))
         return(parameters_ssr_linear_str)
      }
   }
   system.time({param_ssr_matrix = foreach(i = 1:nrow(combinations), .combine = rbind, .packages = c('doParallel', 'minpack.lm')) %dopar% {get_init_param_i(i)}})
   param_ssr_matrix = na.omit(param_ssr_matrix)
   index_selection = which(param_ssr_matrix[,ncol(param_ssr_matrix)] == min(param_ssr_matrix[,ncol(param_ssr_matrix)])) 
   par_str_ini = c(param_ssr_matrix[index_selection,-ncol(param_ssr_matrix)],gamma_grid[combinations[index_selection,1]],quantile(state_variable,probs = c1_grid[combinations[index_selection,2]]),quantile(state_variable,probs = c2_grid[combinations[index_selection,3]]))

   # Estimate by NLS LM
   nls_str_LM = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4
                              +b5*x5+b6*x6 
                              -rho*b1*x7-rho*b2*x8-rho*b3*x9-rho*b4*x10
                              -rho*b5*x11-rho*b6*x12
                              +(1/(1+exp(-(gamma/sd(state_variable[2:obs])^2)*(s_var-c1)*(s_var-c2))))*(b7 + b8*x1 + b9*x2 + b10*x3 + b11*x4 + b12*x5 + b13*x6)
                              -rho*(1/(1+exp(-(gamma/sd(state_variable[2:obs])^2)*(s_var_l-c1)*(s_var_l-c2))))*(b7 + b8*x7 + b9*x8 +b10*x9 + b11*x10 + b12*x11 + b13*x12)
                              , 
                              data = list(y = data_str_inf$inflation[2:obs], y_l = data_str_inf$inflation[1:obs-1], 
                                          x1 = data_str_inf$inflation.l1[2:obs], x2 = data_str_inf$output_gap[2:obs], x3 = data_str_inf$output_gap.l1[2:obs], 
                                          x4 = data_str_inf$real_interest_gap.l2[2:obs], x5 = data_str_inf$rer_gap[2:obs], x6 = data_str_inf$inflation_ext.l2[2:obs],
                                          x7 = data_str_inf$inflation.l1[1:obs-1], x8 = data_str_inf$output_gap[1:obs-1], x9 = data_str_inf$output_gap.l1[1:obs-1], 
                                          x10 = data_str_inf$real_interest_gap.l2[1:obs-1], x11 = data_str_inf$rer_gap[1:obs-1], x12 = data_str_inf$inflation_ext.l2[1:obs-1],
                                          s_var = state_variable[2:obs], s_var_l = state_variable[1:obs-1]) 
                              ,
                              start = list(rho = par_str_ini[1], b0 = par_str_ini[2], b1 = par_str_ini[3], b2 = par_str_ini[4], b3= par_str_ini[5], b4= par_str_ini[6], b5= par_str_ini[7], b6 = par_str_ini[8], 
                                           b7 = par_str_ini[9], b8 = par_str_ini[10], b9 = par_str_ini[11], b10 = par_str_ini[12], b11 = par_str_ini[13], b12 = par_str_ini[14],b13 = par_str_ini[15],
                                           gamma = par_str_ini[16], c1 = par_str_ini[17],c2 = par_str_ini[18]),
                              
                              lower = c(0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,gamma_grid[1],quantile(state_variable,probs = c1_grid[1]),quantile(state_variable,probs = c2_grid[1])),
                              upper = c(Inf, Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf, Inf, quantile(state_variable,probs = c1_grid[length(c_grid)]),quantile(state_variable,probs = c1_grid[length(c2_grid)])),
                              control = nls.lm.control(maxiter = 200), nls.control(minFactor = 1/4096, warnOnly = TRUE), jac = NULL), silent = FALSE)
    
   if(try(summary(nls_str_LM),silent = TRUE)[2] %in% "try-error" || class(nls_str_LM) %in% "try-error"){next} else {
      convergence_message=if(nls_str_LM$convInfo$isConv == TRUE){convergence_message = "converged"} else {"no convergence"}
      #Thing we need to get from the model:
      # Name of the transition variable 
      name_state_variable = names_str_inf[i]
      # Parameters (linear and nonlinear part) 
      parameters_aux = nls_str_LM$m$getPars()
      #Std error of parameter estimates
      std_parameters_aux = summary(nls_str_LM)$coefficients[,2]
      #Transition Function (dated and ordered)
      f_lstr2_aux = 1/(1+exp(-(parameters_aux[length(parameters_aux)-2]/sd(state_variable[2:obs])^2)*(state_variable[2:obs]-parameters_aux[length(parameters_aux)-1])*(state_variable[2:obs]-parameters_aux[length(parameters_aux)])))
         #plot(f_lstr2_aux, type = 'l', col = 'blue')
         #plot(state_variable[2:obs],f_lstr2_aux,  type = 'points')
      f_lstr2_aux_ordered = sort(f_lstr2_aux, decreasing = FALSE)
         #plot(f_lstr2_aux_ordered, col = 'red')
      #Predicted values: 
      predicted_str_aux = predict(nls_str_LM)
      #BIC (to compare against linear models) 
      BIC_str_aux = BIC(nls_str_LM)
      BIC_criteria = if(BIC_str_aux >= BIC_inf_nls){BIC_criteria = "Linear"} else {BIC_criteria = "No Linear"}
      #AIC (to compare against alternative STR especifications)
      AIC_str_aux = AIC(nls_str_LM)
      AIC_Criteria = if(AIC_str_aux >= AIC_inf_nls){AIC_Criteria = "Linear"} else {AIC_Criteria = "No Linear"}
      # Residuals of the model (with 0 at the beginning): 
      residuals_str_aux = resid(nls_str_LM)
      # SSR (necessary for evaluation stage)
      ssr_str_aux = sum((predict(nls_str_LM) - data_str_inf$inflation[2:obs])^2)
      SSR_criteria = if(ssr_str_aux >= ssr_1_inf_nls){SSR_criteria = "Linear"} else {SSR_criteria = "No Linear"}
      #Gradient for autocorrelation test 
      gradient_aux = nls_str_LM$m$gradient()
      #Test of restriction of nonlinear parameters = 0
      F_test_nlparam = ((ssr_1_inf_nls-ssr_str_aux)/(summary(inflation_eq_nls)$df[2] -summary(nls_str_LM)$df[2]))/(ssr_str_aux/(summary(nls_str_LM)$df[2]))
      pval_f_nlparam = pf(F_test_nlparam,summary(inflation_eq_nls)$df[2]-summary(nls_str_LM)$df[2],summary(nls_str_LM)$df[2],lower.tail = FALSE)
      #Test of autocorrelation of residuals 
      #1. Orthogonalize errors: 
      ortho_error_autocorr = lm(residuals_str_aux ~ 1 +gradient_aux)
      reisduals_ortho_str_aux = resid(ortho_error_autocorr)
      #2. Matrix of lagged residuales (vt):
      vt_l_residuals = cbind(reisduals_ortho_str_aux,dplyr::lag(reisduals_ortho_str_aux,1),dplyr::lag(reisduals_ortho_str_aux,2),dplyr::lag(reisduals_ortho_str_aux,3),
                             dplyr::lag(reisduals_ortho_str_aux,4),dplyr::lag(reisduals_ortho_str_aux,5),dplyr::lag(reisduals_ortho_str_aux,6),dplyr::lag(reisduals_ortho_str_aux,7),
                             dplyr::lag(reisduals_ortho_str_aux,9),dplyr::lag(reisduals_ortho_str_aux,10),dplyr::lag(reisduals_ortho_str_aux,11),dplyr::lag(reisduals_ortho_str_aux,12))
      vt_l_residuals[is.na(vt_l_residuals)] <- 0
      vt_l_residuals = vt_l_residuals[,c(-1)]
      #3. Regression result and LM tests (F and Chi) 
      auto_corr_test = lm(residuals_str_aux ~ gradient_aux + vt_l_residuals)
      ssr_alt = sum((predict(auto_corr_test)-residuals_str_aux)^2) 
      F_corr_str = ((ssr_str_aux - ssr_alt)/(summary(nls_str_LM)$df[2]-auto_corr_test$df.residual))/(ssr_alt/(auto_corr_test$df.residual))
      pval_f_corr = pf(F_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual,auto_corr_test$df.residual,lower.tail = FALSE)
      chi_corr_str = length(residuals_str_aux)*((ssr_str_aux - ssr_alt)/ssr_str_aux) 
      pval_chi_corr = pchisq(chi_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual, lower.tail = FALSE) 
      #Result of Jarque-Bera and Shapiro Test 
      jb_test_aux = jarque.bera.test(residuals_str_aux)$p.value
      shapiro_aux = shapiro.test(residuals_str_aux)$p.value
      #Ramsey Test
      ols_str_ramsey = lm(residuals_str_aux ~ gradient_aux + predicted_str_aux + predicted_str_aux^2 + predicted_str_aux^3)
      ssr_ols_str_ramsey = sum((predict(ols_str_ramsey) - residuals_str_aux)^2)
      #F test statistic:
      f_ramsey_ols = ((ssr_str_aux - ssr_ols_str_ramsey)/(summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual))/(ssr_ols_str_ramsey/(ols_str_ramsey$df.residual))
      pval_f_str_ramsey = pf(f_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,ols_str_ramsey$df.residual, lower.tail = FALSE)
      #Chi test statistic 
      chi_ramsey_ols = length(residuals_str_aux)*(ssr_str_aux - ssr_ols_str_ramsey)/ssr_str_aux
      pval_chi_str_ramsey = pchisq(chi_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,lower.tail = FALSE)
      #No additive Linearity Test - LM: 
      no_ad_nolin_results_direct = NULL
      for(k in eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct))){
         state_2 =  data_str_inf[2:nrow(data_str_inf),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux + eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_direct = rbind(no_ad_nolin_results_direct,no_ad_nolin_results_aux)
      }
      
      no_ad_nolin_results_state = NULL
      for(k in eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state))){
         state_2 =  data_str_inf[2:nrow(data_str_inf),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux +state_2+eval(state_2^2)+eval(state_2^3)+ eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_state = rbind(no_ad_nolin_results_state,no_ad_nolin_results_aux)
      }
      no_ad_nolin_results_direct = na.omit(no_ad_nolin_results_direct)
      no_ad_nolin_results_state = na.omit(no_ad_nolin_results_state)
      if(any(no_ad_nolin_results_direct[,2] <= 0.05) == TRUE){no_ad_d_nl_5_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_direct[,2] <= 0.01) == TRUE){no_ad_d_nl_1_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_1_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.05) == TRUE){no_ad_s_nl_5_pc = length(which(no_ad_nolin_results_state[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.01) == TRUE){no_ad_s_nl_1_pc = length(which(no_ad_nolin_results_state[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_1_pc="No Additive Nonlinearity"}
      #Parameter Constancy Test - LM: cannot be performed do to small sample size problem (low DF)
      #Arch test: 
      noarch_test =  ArchTest(residuals_str_aux, lags = 12)
      arch_chi_statistic = noarch_test$statistic
      arch_chi_pval = noarch_test$p.value
      
      #Final Output
      model_info = rbind(convergence_message,BIC_str_aux,BIC_criteria, AIC_str_aux, AIC_Criteria, ssr_str_aux,SSR_criteria,F_corr_str,pval_f_corr,chi_corr_str,pval_chi_corr,
                         jb_test_aux,shapiro_aux,f_ramsey_ols,pval_f_str_ramsey,chi_ramsey_ols,pval_chi_str_ramsey,no_ad_d_nl_5_pc,no_ad_d_nl_1_pc,no_ad_s_nl_5_pc,no_ad_s_nl_1_pc,arch_chi_statistic,arch_chi_pval, 
                         F_test_nlparam, pval_f_nlparam) 
      colnames(model_info) = name_state_variable
      rownames(model_info) = c("Convergence", "BIC", "BIC Criteria", "AIC", "AIC Criteria", "SSR", "SSR Criteria", "F test corr", "Pval F test", "Chi2 test corr", "Pval Chi2 test",
                               "JB test", "Shapiro test", "F test Ramsey", "Pval F Ramsey", "Chi2 test Ramsey", "Pval Chi2 Ramsey", "No Ad NL D 5%","No Ad D NL 1%","No Ad NL S 5%","No Ad S NL 1%", "ARCH test stat", "ARCH test pval", 
                               "F test nl", "Pval F nl")
      
      Transition_function = cbind(f_lstr1_aux, f_lstr1_aux_ordered, state_variable[2:obs])
      colnames(Transition_function) = c("LSTR2", "LSTR2 sorted", "State_Variable")
      final_output = list(Model_Info = model_info,Model = nls_str_LM, Transition_function = Transition_function)
      LSTR2_INF[[eval(i-ncol(data_str_inf_model))]] <- final_output
   }
}
names(LSTR2_INF) <- names_str_inf[eval(ncol(data_str_inf_model)+1):length(names_str_inf)]
proc.time() - ptm 

#Results: of model info 
LSTR2_INF_model_results = NULL
for(j in 1:length(LSTR2_INF)){aux_info = LSTR2_INF[[j]]$Model_Info
LSTR2_INF_model_results = cbind(LSTR2_INF_model_results,aux_info)
}
#Selected models According to BIC/AIC criteria,  
LSTR2_INF_selected_models = LSTR2_INF_model_results[,
                                                  which(LSTR2_INF_model_results['Convergence',] == "converged" & 
                                                           (LSTR2_INF_model_results['BIC Criteria',] == "No Linear" |  LSTR2_INF_model_results['AIC Criteria',] == "No Linear")
                                                        & (LSTR2_INF_model_results['JB test',] >= 0.05 | LSTR2_INF_model_results['Shapiro test',] >= 0.05)
                                                        &  LSTR2_INF_model_results['Pval F nl',] <= 0.05 &  LSTR2_INF_model_results['Pval F test',] >= 0.05 &  LSTR2_INF_model_results['Pval Chi2 test',] >= 0.05
                                                        &  LSTR2_INF_model_results['Pval F Ramsey',] >= 0.05 & LSTR2_INF_model_results['Pval Chi2 Ramsey',] >= 0.05
                                                        & LSTR2_INF_model_results['ARCH test pval',] >= 0.05 & (LSTR2_INF_model_results['No Ad D NL 1%',] == "No Additive Nonlinearity" | LSTR2_INF_model_results['No Ad D NL 1%',] < 0.25) 
                                                        & (LSTR2_INF_model_results['No Ad S NL 1%',] == "No Additive Nonlinearity" | LSTR2_INF_model_results['No Ad S NL 1%',] < 0.25))]

write.csv(LSTR2_INF_selected_models, "LSTR2_INF_selected.csv")

# ------- Output Gap Equation ------------- #
# Selection of variables where LM test detected non-linearity
names_outputgap_direct = colnames(LM_output_gap_direct)
names_outputgap_state = colnames(LM_outputgap_state)
selected_variables_direct = names_outputgap_direct[which(LM_output_gap_direct[c("LM result"),] == "Nonlinear")]
selected_variables_state = names_outputgap_state[which(LM_outputgap_state[c("LM result"),] == "Nonlinear")]
#selected_variables_state = names_outputgap_state[which(LM_outputgap_state[c("LM result"),] == "Nonlinear" && LM_outputgap_state[c("STR model F"),] == "LSTR2-ESTR")]
selected_variables = c(selected_variables_direct,selected_variables_state)

data_str_outputgap_state = data_state_outputgap[,selected_variables] 
data_str_outputgap_state = xts(data_str_outputgap_state, index_dates_data_state_final)
data_str_outputgap_model = as.data.frame(cbind(data_state_final$output_gap, data_state_final$resid_outputgap, data_state_final$output_gap.l1, data_state_final$real_interest_gap.l2,data_state_final$rer_gap, data_state_final$us_output_gap,
                                         data_state_final$brent_gap))
data_str_outputgap_model = xts(data_str_outputgap_model, index_dates_data_state_final)
data_str_outputgap = cbind.xts(data_str_outputgap_model, data_str_outputgap_state)

obs = nrow(data_str_outputgap)
data_str_outputgap = as.data.frame(data_str_outputgap)
names_str_outputgap = colnames(data_str_outputgap)
ncol(data_str_outputgap[,eval(ncol(data_str_outputgap_model)+1):eval(ncol(data_str_outputgap_model)+ncol(data_str_outputgap_state))])

#Steps (for each nonlinear or state variable)
#1. Calculate Grid for Gamma and C parameters
#2. USe the Grid to compute the transition function (logistic/exponential)
#3. Add transition function to the model and estimate model by OLS - NLS
#4. Obtain SSR for each iteration, and choose the parameters that minimize it
#5. Use parameters in 4 as starting values to NLS estimate of STR


LSTR1_OUTPUTGAP = list()
ptm <- proc.time()
#For Loop for estimates of LSTR1 function 
for(i in eval(ncol(data_str_outputgap_model)+1):eval(ncol(data_str_outputgap_model)+ncol(data_str_outputgap_state))){
   #Declare state variable
   state_variable = data_str_outputgap[,i]  
   #Step 1: Grid for Gamma and C parameters
   gamma_grid = seq(0.1*sd(state_variable[2:obs]), 32*sd(state_variable[2:obs]), length.out = 32) #Reescaled by transition variable sd
   c_grid = seq(0.1, 0.95, 0.02) #From 10th percentile to 90th percentile
   combinations = as.data.frame(crossing(gamma = 1:length(gamma_grid), c = 1:length(c_grid)))
   dim(combinations)
   
   #LSTR 1 estimation
   ssr_linear_str = NULL
   parameters_linear_str = NULL
   ptm <- proc.time()
   get_init_param_o <- function(x){
      gamma_cond <- gamma_grid[combinations[x,1]]
      c_cond <-quantile(state_variable, probs = c_grid[combinations[x,2]]) #State variable transition point
      # Compute transition function for the state variable
      f_lstr1 = 1/(1+exp(-(gamma_cond/sd(state_variable[2:obs]))*(state_variable-c_cond))) #This becomes a variable 
      # Estimate model by OLS (without rho term) to get initial parameters for NLS estimation
      ols_aux = lm(output_gap  ~ 1 + output_gap.l1 + real_interest_gap.l2 + rer_gap + us_output_gap + brent_gap             
                   + f_lstr1 + eval(output_gap.l1*f_lstr1) + eval(real_interest_gap.l2*f_lstr1) + eval(rer_gap*f_lstr1) + eval(us_output_gap*f_lstr1) + eval(brent_gap*f_lstr1),
                   data = cbind(data_str_outputgap,f_lstr1))
      par_ols = coef(ols_aux)
      par_ols[is.na(par_ols)] <- 0
      # Add variable to the estimation of linear model (by NLS):
      #From x11 to x15 are the STR interactions of model variables with state variables
      #From x16 to x20 are lagged interactions with rho term
      
      linear_str = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4 +b5*x5
                             -rho*b1*x6-rho*b2*x7-rho*b3*x8-rho*b4*x9-rho*b5*x10
                             +b6*x11 + b7*x12 + b8*x13 + b9*x14 + b10*x15 + b11*x16
                             - rho*b6*x17 - rho*b7*x18 - rho*b8*x19- rho*b9*x20 - rho*b10*x21 - rho*b11*x22
                             , 
                             start = list(rho = 0.03, b0 = -0.0015, b1 = 0.88, b2 = -0.023, b3= -0.084, b4= 0.23, b5= 0.004, 
                                          b6 = par_ols[7], b7 = par_ols[8], b8 = par_ols[9], b9 = par_ols[10], b10 = par_ols[11], b11 = par_ols[12]),
                             data = list(y = data_str_outputgap$output_gap[2:obs], y_l = data_str_outputgap$output_gap[1:obs-1], x1 = data_str_outputgap$output_gap.l1[2:obs],
                                         x2 = data_str_outputgap$real_interest_gap.l2[2:obs], x3 = data_str_outputgap$us_output_gap[2:obs], x4 = data_str_outputgap$rer_gap[2:obs], x5 = data_str_outputgap$brent_gap[2:obs], 
                                         x6 = data_str_outputgap$output_gap.l1[1:obs-1], x7 = data_str_outputgap$real_interest_gap.l2[1:obs-1], x8 = data_str_outputgap$us_output_gap[1:obs-1], x9 = data_str_outputgap$rer_gap[1:obs-1], x10 = data_str_outputgap$brent_gap[1:obs-1],
                                         x11 = f_lstr1[2:obs], x12 = data_str_outputgap$output_gap.l1[2:obs]*f_lstr1[2:obs], x13 = data_str_outputgap$real_interest_gap.l2[2:obs]*f_lstr1[2:obs], x14 = data_str_outputgap$us_output_gap[2:obs]*f_lstr1[2:obs], x15 = data_str_outputgap$rer_gap[2:obs]*f_lstr1[2:obs], x16 = data_str_outputgap$brent_gap[2:obs]*f_lstr1[2:obs], 
                                         x17 = f_lstr1[1:obs-1], x18 = data_str_outputgap$output_gap.l1[1:obs-1]*f_lstr1[1:obs-1], x19 = data_str_outputgap$real_interest_gap.l2[1:obs-1]*f_lstr1[1:obs-1], x20 = data_str_outputgap$us_output_gap[1:obs-1]*f_lstr1[1:obs-1], x21 = data_str_outputgap$rer_gap[1:obs-1]*f_lstr1[1:obs-1], x22 = data_str_outputgap$brent_gap[1:obs-1]*f_lstr1[1:obs-1]
                                         
                             ), 
                             trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200) , nls.control(minFactor = 1/4096, warnOnly = TRUE)), silent = FALSE)
      if(try(summary(linear_str),silent = TRUE)[2] %in% "try-error" || class(linear_str) %in% "try-error"){return(c(rep(NA, 14)))} else {
         # Get the Residual Sum of Squares of the Regression
         ssr_aux = sum((predict(linear_str) - data_str_outputgap$output_gap[2:obs]) ^ 2)  ## residual sum of squares
         #ssr_linear_str = rbind(ssr_linear_str, ssr_aux)
         # Get the model parameters
         parameters_aux = linear_str$m$getPars()
         parameters_ssr_linear_str = rbind(c(), c(parameters_aux, ssr_aux))
         return(parameters_ssr_linear_str)
      }
   }
   system.time({param_ssr_matrix = foreach(i = 1:nrow(combinations), .combine = rbind, .packages = c('doParallel', 'minpack.lm')) %dopar% {get_init_param_o(i)}})
   param_ssr_matrix = na.omit(param_ssr_matrix)
   index_selection = which(param_ssr_matrix[,ncol(param_ssr_matrix)] == min(param_ssr_matrix[,ncol(param_ssr_matrix)])) 
   par_str_ini = c(param_ssr_matrix[index_selection,-ncol(param_ssr_matrix)],gamma_grid[combinations[index_selection,1]],quantile(state_variable,probs = c_grid[combinations[index_selection,2]]))

   # Estimate by NLS LM
   nls_str_LM = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4+b5*x5
                              -rho*b1*x6-rho*b2*x7-rho*b3*x8-rho*b4*x9-rho*b5*x10
                              +(1/(1+exp(-(gamma/sd(state_variable[2:obs]))*(s_var-c))))*(b6+b7*x1+b8*x2+b9*x3+b10*x4+b11*x5)
                              -rho*(1/(1+exp(-(gamma/sd(state_variable[2:obs]))*(s_var_l-c))))*(b6 + b7*x6+b8*x7+b9*x8+b10*x9+b11*x10)
                              , 
                          data = list(y = data_str_outputgap$output_gap[2:obs], y_l = data_str_outputgap$output_gap[1:obs-1], x1 = data_str_outputgap$output_gap.l1[2:obs],
                                      x2 = data_str_outputgap$real_interest_gap.l2[2:obs], x3 = data_str_outputgap$us_output_gap[2:obs], x4 = data_str_outputgap$rer_gap[2:obs], x5 = data_str_outputgap$brent_gap[2:obs], 
                                      x6 = data_str_outputgap$output_gap.l1[1:obs-1], x7 = data_str_outputgap$real_interest_gap.l2[1:obs-1], x8 = data_str_outputgap$us_output_gap[1:obs-1], x9 = data_str_outputgap$rer_gap[1:obs-1], x10 = data_str_outputgap$brent_gap[1:obs-1],
                                      s_var = state_variable[2:obs], s_var_l = state_variable[1:obs-1]) 
                              ,
                              start = list(rho = par_str_ini[1], b0 = par_str_ini[2], b1 = par_str_ini[3], b2 = par_str_ini[4], b3= par_str_ini[5], b4= par_str_ini[6], b5= par_str_ini[7], b6 = par_str_ini[8], 
                                           b7 = par_str_ini[9], b8 = par_str_ini[10], b9 = par_str_ini[11], b10 = par_str_ini[12], b11 = par_str_ini[13], gamma = par_str_ini[14], c = par_str_ini[15]),
                              
                              lower = c(0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf, -Inf,-Inf,gamma_grid[1],quantile(state_variable,probs = c_grid[1])),
                              upper = c(Inf, Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,quantile(state_variable,probs = c_grid[length(c_grid)])),
                              control = nls.lm.control(maxiter = 200), nls.control(minFactor = 1/4096, warnOnly = TRUE), jac = NULL), silent = TRUE)
   
   if(try(summary(nls_str_LM),silent = TRUE)[2] %in% "try-error" || class(nls_str_LM) %in% "try-error"){next} else { 
      convergence_message=if(nls_str_LM$convInfo$isConv == TRUE){convergence_message = "converged"} else {"no convergence"}
      #Thing we need to get from the model:
      # Name of the transition variable 
      name_state_variable = names_str_outputgap[i]
      # Parameters (linear and nonlinear part) 
      parameters_aux = nls_str_LM$m$getPars()
      #Std error of parameter estimates
      std_parameters_aux = summary(nls_str_LM)$coefficients[,2]
      #Transition Function (dated and ordered)
      f_lstr1_aux = 1/(1+exp(-(parameters_aux[length(parameters_aux)-1]/sd(state_variable[2:obs]))*(state_variable[2:obs]-parameters_aux[length(parameters_aux)])))
      #plot(f_lstr1_aux, type = 'l', col = 'blue')
      #plot(state_variable[2:obs],f_lstr1_aux)
      f_lstr1_aux_ordered = sort(f_lstr1_aux, decreasing = FALSE)
      #plot(f_lstr1_aux_ordered, col = 'red')
      #Predicted values: 
      predicted_str_aux = predict(nls_str_LM)
      #BIC (to compare against linear models) 
      BIC_str_aux = BIC(nls_str_LM)
      BIC_criteria = if(BIC_str_aux >= BIC_outputgap_nls){BIC_criteria = "Linear"} else {BIC_criteria = "No Linear"}
      #AIC (to compare against alternative STR especifications)
      AIC_str_aux = AIC(nls_str_LM)
      AIC_Criteria = if(AIC_str_aux >= AIC_outputgap_nls){AIC_Criteria = "Linear"} else {AIC_Criteria = "No Linear"}
      # Residuals of the model (with 0 at the beginning): 
      residuals_str_aux = resid(nls_str_LM)
      # SSR (necessary for evaluation stage)
      ssr_str_aux = sum((predict(nls_str_LM) - data_str_outputgap$output_gap[2:obs])^2)
      SSR_criteria = if(ssr_str_aux >= ssr_1_outputgap_nls){SSR_criteria = "Linear"} else {SSR_criteria = "No Linear"}
      #Gradient for autocorrelation test 
      gradient_aux = nls_str_LM$m$gradient()
      #Test of restriction of nonlinear parameters = 0
      F_test_nlparam = ((ssr_1_outputgap_nls-ssr_str_aux)/(summary(outputgap_eq_nls)$df[2] -summary(nls_str_LM)$df[2]))/(ssr_str_aux/(summary(nls_str_LM)$df[2]))
      pval_f_nlparam = pf(F_test_nlparam,summary(outputgap_eq_nls)$df[2]-summary(nls_str_LM)$df[2],summary(nls_str_LM)$df[2],lower.tail = FALSE)
      #Test of autocorrelation of residuals 
      #1. Orthogonalize errors: 
      ortho_error_autocorr = lm(residuals_str_aux ~ 1 +gradient_aux)
      reisduals_ortho_str_aux = resid(ortho_error_autocorr)
      #2. Matrix of lagged residuales (vt):
      vt_l_residuals = cbind(reisduals_ortho_str_aux,dplyr::lag(reisduals_ortho_str_aux,1),dplyr::lag(reisduals_ortho_str_aux,2),dplyr::lag(reisduals_ortho_str_aux,3),
                             dplyr::lag(reisduals_ortho_str_aux,4),dplyr::lag(reisduals_ortho_str_aux,5),dplyr::lag(reisduals_ortho_str_aux,6),dplyr::lag(reisduals_ortho_str_aux,7),
                             dplyr::lag(reisduals_ortho_str_aux,9),dplyr::lag(reisduals_ortho_str_aux,10),dplyr::lag(reisduals_ortho_str_aux,11),dplyr::lag(reisduals_ortho_str_aux,12))
      vt_l_residuals[is.na(vt_l_residuals)] <- 0
      vt_l_residuals = vt_l_residuals[,c(-1)]
      #3. Regression result and LM tests (F and Chi) 
      auto_corr_test = lm(residuals_str_aux ~ gradient_aux + vt_l_residuals)
      ssr_alt = sum((predict(auto_corr_test)-residuals_str_aux)^2) 
      F_corr_str = ((ssr_str_aux - ssr_alt)/(summary(nls_str_LM)$df[2]-auto_corr_test$df.residual))/(ssr_alt/(auto_corr_test$df.residual))
      pval_f_corr = pf(F_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual,auto_corr_test$df.residual,lower.tail = FALSE)
      chi_corr_str = length(residuals_str_aux)*((ssr_str_aux - ssr_alt)/ssr_str_aux) 
      pval_chi_corr = pchisq(chi_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual, lower.tail = FALSE) 
      #Result of Jarque-Bera and Shapiro Test 
      jb_test_aux = jarque.bera.test(residuals_str_aux)$p.value
      shapiro_aux = shapiro.test(residuals_str_aux)$p.value
      #Ramsey Test
      ols_str_ramsey = lm(residuals_str_aux ~ gradient_aux + predicted_str_aux + predicted_str_aux^2 + predicted_str_aux^3)
      ssr_ols_str_ramsey = sum((predict(ols_str_ramsey) - residuals_str_aux)^2)
      #F test statistic:
      f_ramsey_ols = ((ssr_str_aux - ssr_ols_str_ramsey)/(summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual))/(ssr_ols_str_ramsey/(ols_str_ramsey$df.residual))
      pval_f_str_ramsey = pf(f_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,ols_str_ramsey$df.residual, lower.tail = FALSE)
      #Chi test statistic 
      chi_ramsey_ols = length(residuals_str_aux)*((ssr_str_aux - ssr_ols_str_ramsey)/ssr_str_aux)
      pval_chi_str_ramsey = pchisq(chi_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,lower.tail = FALSE)
      #No additive Linearity Test - LM: 
      no_ad_nolin_results_direct = NULL
      for(k in eval(ncol(data_str_outputgap_model)+1):eval(ncol(data_str_outputgap_model)+length(selected_variables_direct))){
         state_2 =  data_str_outputgap[2:nrow(data_str_outputgap),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux + eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_direct = rbind(no_ad_nolin_results_direct,no_ad_nolin_results_aux)
      }
      
      no_ad_nolin_results_state = NULL
      for(k in eval(ncol(data_str_outputgap_model)+length(selected_variables_direct)+1):eval(ncol(data_str_outputgap_model)+ncol(data_str_outputgap_state))){
         state_2 =  data_str_outputgap[2:nrow(data_str_outputgap),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux +state_2+eval(state_2^2)+eval(state_2^3)+ eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_state = rbind(no_ad_nolin_results_state,no_ad_nolin_results_aux)
      }
      no_ad_nolin_results_direct = na.omit(no_ad_nolin_results_direct)
      no_ad_nolin_results_state = na.omit(no_ad_nolin_results_state)
      if(any(no_ad_nolin_results_direct[,2] <= 0.05) == TRUE){no_ad_d_nl_5_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_direct[,2] <= 0.01) == TRUE){no_ad_d_nl_1_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_1_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.05) == TRUE){no_ad_s_nl_5_pc = length(which(no_ad_nolin_results_state[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.01) == TRUE){no_ad_s_nl_1_pc = length(which(no_ad_nolin_results_state[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_1_pc="No Additive Nonlinearity"}
      #Parameter Constancy Test - LM: cannot be performed do to small sample size problem (low DF)
      #Arch test: 
      noarch_test =  ArchTest(residuals_str_aux, lags = 12)
      arch_chi_statistic = noarch_test$statistic
      arch_chi_pval = noarch_test$p.value
      
      #Final Output
      model_info = rbind(convergence_message,BIC_str_aux,BIC_criteria, AIC_str_aux, AIC_Criteria, ssr_str_aux,SSR_criteria,F_corr_str,pval_f_corr,chi_corr_str,pval_chi_corr,
                         jb_test_aux,shapiro_aux,f_ramsey_ols,pval_f_str_ramsey,chi_ramsey_ols,pval_chi_str_ramsey,no_ad_d_nl_5_pc,no_ad_d_nl_1_pc,no_ad_s_nl_5_pc,no_ad_s_nl_1_pc,arch_chi_statistic,arch_chi_pval, 
                         F_test_nlparam, pval_f_nlparam) 
      colnames(model_info) = name_state_variable
      rownames(model_info) = c("Convergence", "BIC", "BIC Criteria", "AIC", "AIC Criteria", "SSR", "SSR Criteria", "F test corr", "Pval F test", "Chi2 test corr", "Pval Chi2 test",
                               "JB test", "Shapiro test", "F test Ramsey", "Pval F Ramsey", "Chi2 test Ramsey", "Pval Chi2 Ramsey", "No Ad NL D 5%","No Ad D NL 1%","No Ad NL S 5%","No Ad S NL 1%", "ARCH test stat", "ARCH test pval", 
                               "F param nl", "Pval F param nl")
      Transition_function = cbind(f_lstr1_aux, f_lstr1_aux_ordered, state_variable[2:obs])
      colnames(Transition_function) = c("LSTR1", "LSTR1 sorted", "State_Variable")
      final_output = list(Model_Info = model_info,Model = nls_str_LM, Transition_function = Transition_function)
      LSTR1_OUTPUTGAP[[eval(i-ncol(data_str_outputgap_model))]] <- final_output
   }
}
names(LSTR1_OUTPUTGAP) <- names_str_outputgap[eval(ncol(data_str_outputgap_model)+1):eval(length(names_str_outputgap)-1)]
proc.time() - ptm 

aux_lstr1_outpugap_2 = NULL
for(i in 1:70){
   if(is.null(LSTR1_OUTPUTGAP[[i]]$Transition_function[1,3])){aux = NA} else {
   aux = LSTR1_OUTPUTGAP[[i]]$Transition_function[1,3]}
   aux_lstr1_outpugap_2  = rbind(aux_lstr1_outpugap_2,aux)
}
aux_lstr1_outpugap_2 = c(aux_lstr1_outpugap_2, NA)
diff_aux = as.data.frame(aux_lstr1_outpugap - aux_lstr1_outpugap_2)

#Results: of model info 
LSTR1_OUTPUTGAP_model_results = NULL
for(j in 1:length(LSTR1_OUTPUTGAP)){aux_info = LSTR1_OUTPUTGAP[[j]]$Model_Info
LSTR1_OUTPUTGAP_model_results = cbind(LSTR1_OUTPUTGAP_model_results,aux_info)
}
#Selected models According to BIC/AIC criteria,  
LSTR1_OUTPUTGAP_selected_models = LSTR1_OUTPUTGAP_model_results[,
                                                    which(LSTR1_OUTPUTGAP_model_results['Convergence',] == "converged" & 
                                                             (LSTR1_OUTPUTGAP_model_results['BIC Criteria',] == "No Linear" |  LSTR1_OUTPUTGAP_model_results['AIC Criteria',] == "No Linear")
                                                          & (LSTR1_OUTPUTGAP_model_results['JB test',] >= 0.05 | LSTR1_OUTPUTGAP_model_results['Shapiro test',] >= 0.05)
                                                          &  LSTR1_OUTPUTGAP_model_results['Pval F param nl',] <= 0.05 &  LSTR1_OUTPUTGAP_model_results['Pval F test',] >= 0.05 &  LSTR1_OUTPUTGAP_model_results['Pval Chi2 test',] >= 0.05
                                                          &  LSTR1_OUTPUTGAP_model_results['Pval F Ramsey',] >= 0.05 & LSTR1_OUTPUTGAP_model_results['Pval Chi2 Ramsey',] >= 0.05
                                                          & LSTR1_OUTPUTGAP_model_results['ARCH test pval',] >= 0.05 & (LSTR1_OUTPUTGAP_model_results['No Ad D NL 1%',] == "No Additive Nonlinearity" | LSTR1_OUTPUTGAP_model_results['No Ad D NL 1%',] < 0.25) 
                                                          & (LSTR1_OUTPUTGAP_model_results['No Ad S NL 1%',] == "No Additive Nonlinearity" | LSTR1_OUTPUTGAP_model_results['No Ad S NL 1%',] < 0.25))]

write.csv(LSTR1_OUTPUTGAP_selected_models, "LSTR1_OUTPUTGAP_selected.csv")

# Selection of variables where LM test detected non-linearity and the Wald test suggested ESTR or LSTR2 model specification: 
names_outputgap_direct = colnames(LM_output_gap_direct)
names_outputgap_state = colnames(LM_outputgap_state)
selected_variables_direct = names_outputgap_direct[which(LM_output_gap_direct[c("LM result"),] == "Nonlinear" & LM_output_gap_direct[c("STR model F"),] == "LSTR2-ESTR")]
selected_variables_state = names_outputgap_state[which(LM_outputgap_state[c("LM result"),] == "Nonlinear" & LM_outputgap_state[c("STR model F"),] == "LSTR2-ESTR")]
selected_variables = c(selected_variables_direct,selected_variables_state)

data_str_outputgap_state = data_state_outputgap[,selected_variables] 
data_str_outputgap_state = xts(data_str_outputgap_state, index_dates_data_state_final)
data_str_outputgap_model = as.data.frame(cbind(data_state_final$output_gap, data_state_final$resid_outputgap, data_state_final$output_gap.l1, data_state_final$real_interest_gap.l2,data_state_final$rer_gap, data_state_final$us_output_gap,
                                               data_state_final$brent_gap))
data_str_outputgap_model = xts(data_str_outputgap_model, index_dates_data_state_final)
data_str_outputgap = cbind.xts(data_str_outputgap_model, data_str_outputgap_state)

obs = nrow(data_str_outputgap)
data_str_outputgap = as.data.frame(data_str_outputgap)
names_str_outputgap = colnames(data_str_outputgap)

ESTR_OUTPUTGAP = list()
ptm <- proc.time()
#For Loop for estimates of LSTR1 function 
for(i in eval(ncol(data_str_outputgap_model)+1):eval(ncol(data_str_outputgap_model)+ncol(data_str_outputgap_state))){
   #Declare state variable
   state_variable = data_str_outputgap[,i]  
   #Step 1: Grid for Gamma and C parameters
   gamma_grid = seq(0.1*sd(state_variable[2:obs]), 32*sd(state_variable[2:obs]), length.out = 32) #Reescaled by transition variable sd
   c_grid = seq(0.1, 0.95, 0.02) #From 10th percentile to 90th percentile
   combinations = as.data.frame(crossing(gamma = 1:length(gamma_grid), c = 1:length(c_grid)))
   dim(combinations)
   
   ssr_linear_str = NULL
   parameters_linear_str = NULL
   #ESTR  estimation
   get_init_param_o <- function(x){
      #Combinations of parameters
      gamma_cond <- gamma_grid[combinations[x,1]]
      c_cond <-quantile(state_variable, probs = c_grid[combinations[x,2]]) #State variable transition point
      # Compute transition function for the state variable
      f_estr1 = 1-exp(-(gamma_cond/sd(state_variable[2:obs]))*((state_variable-c_cond)^2)) #This becomes a variable 
      # Estimate model by OLS (without rho term) to get initial parameters for NLS estimation
      ols_aux = lm(output_gap  ~ 1 + output_gap.l1 + real_interest_gap.l2 + rer_gap + us_output_gap + brent_gap             
                   +f_estr1+eval(output_gap.l1*f_estr1) + eval(real_interest_gap.l2*f_estr1) + eval(rer_gap*f_estr1) + eval(us_output_gap*f_estr1) + eval(brent_gap*f_estr1),
                   data = data_str_outputgap)
      par_ols = coef(ols_aux)
      par_ols[is.na(par_ols)]<- 0
      
      # Add variable to the estimation of linear model (by NLS):
      #From x11 to x15 are the STR interactions of model variables with state variables
      #From x16 to x20 are lagged interactions with rho term
      
      linear_str = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4 +b5*x5
                             -rho*b1*x6-rho*b2*x7-rho*b3*x8-rho*b4*x9-rho*b5*x10
                             +b6*x11 + b7*x12 + b8*x13 + b9*x14 + b10*x15 +b11*x16
                             - rho*b6*x17 - rho*b7*x18 - rho*b8*x19- rho*b9*x20 - rho*b10*x21 - rho*b11*x22
                             , 
                             start = list(rho = 0.03, b0 = -0.0015, b1 = 0.88, b2 = -0.023, b3= -0.084, b4= 0.23, b5= 0.004, 
                                          b6 = par_ols[7], b7 = par_ols[8], b8 = par_ols[9], b9 = par_ols[10], b10 = par_ols[11], b11 = par_ols[12]),
                             data = list(y = data_str_outputgap$output_gap[2:obs], y_l = data_str_outputgap$output_gap[1:obs-1], x1 = data_str_outputgap$output_gap.l1[2:obs],
                                         x2 = data_str_outputgap$real_interest_gap.l2[2:obs], x3 = data_str_outputgap$us_output_gap[2:obs], x4 = data_str_outputgap$rer_gap[2:obs], x5 = data_str_outputgap$brent_gap[2:obs], 
                                         x6 = data_str_outputgap$output_gap.l1[1:obs-1], x7 = data_str_outputgap$real_interest_gap.l2[1:obs-1], x8 = data_str_outputgap$us_output_gap[1:obs-1], x9 = data_str_outputgap$rer_gap[1:obs-1], x10 = data_str_outputgap$brent_gap[1:obs-1],
                                         x11 = f_estr1[2:obs], x12  = data_str_outputgap$output_gap.l1[2:obs]*f_estr1[2:obs], x13 = data_str_outputgap$real_interest_gap.l2[2:obs]*f_estr1[2:obs], x14 = data_str_outputgap$us_output_gap[2:obs]*f_estr1[2:obs], x15 = data_str_outputgap$rer_gap[2:obs]*f_estr1[2:obs], x16 = data_str_outputgap$brent_gap[2:obs]*f_estr1[2:obs], 
                                         x17 = f_estr1[1:obs-1],x18 = data_str_outputgap$output_gap.l1[1:obs-1]*f_estr1[1:obs-1], x19 = data_str_outputgap$real_interest_gap.l2[1:obs-1]*f_estr1[1:obs-1], x20 = data_str_outputgap$us_output_gap[1:obs-1]*f_estr1[1:obs-1], x21 = data_str_outputgap$rer_gap[1:obs-1]*f_estr1[1:obs-1], x22 = data_str_outputgap$brent_gap[1:obs-1]*f_estr1[1:obs-1]
                                         
                             ), 
                             trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200) , nls.control(minFactor = 1/4096, warnOnly = TRUE)), silent = FALSE)
      if(try(summary(linear_str),silent = TRUE)[2] %in% "try-error" || class(linear_str) %in% "try-error"){return(c(rep(NA, 14)))} else {
         # Get Sum of Square Residuals
         ssr_aux = sum((predict(linear_str) - data_str_outputgap$output_gap[2:obs]) ^ 2)  ## residual sum of squares
         # Get the model parameters
         parameters_aux = linear_str$m$getPars()
         parameters_ssr_linear_str = rbind(c(), c(parameters_aux, ssr_aux))
         return(parameters_ssr_linear_str)
      }
   }
   system.time({param_ssr_matrix = foreach(i = 1:nrow(combinations), .combine = rbind, .packages = c('doParallel', 'minpack.lm')) %dopar% {get_init_param_o(i)}})
   param_ssr_matrix = na.omit(param_ssr_matrix)
   index_selection = which(param_ssr_matrix[,ncol(param_ssr_matrix)] == min(param_ssr_matrix[,ncol(param_ssr_matrix)])) 
   par_str_ini = c(param_ssr_matrix[index_selection,-ncol(param_ssr_matrix)],gamma_grid[combinations[index_selection,1]],quantile(state_variable,probs = c_grid[combinations[index_selection,2]]))
   
   # Estimate by NLS LM
   nls_str_LM = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4+b5*x5
                          -rho*b1*x6-rho*b2*x7-rho*b3*x8-rho*b4*x9-rho*b5*x10
                          +(1-exp(-(gamma/sd(state_variable[2:obs]))*((s_var-c)^2)))*(b6 +b7*x1+b8*x2+b9*x3+b10*x4+b11*x5)
                          -rho*(1-exp(-(gamma/sd(state_variable[2:obs]))*((s_var_l-c)^2)))*(b6 +b7*x6+b8*x7+b9*x8+b10*x9+b11*x10)
                          , 
                          data = list(y = data_str_outputgap$output_gap[2:obs], y_l = data_str_outputgap$output_gap[1:obs-1], x1 = data_str_outputgap$output_gap.l1[2:obs],
                                      x2 = data_str_outputgap$real_interest_gap.l2[2:obs], x3 = data_str_outputgap$us_output_gap[2:obs], x4 = data_str_outputgap$rer_gap[2:obs], x5 = data_str_outputgap$brent_gap[2:obs], 
                                      x6 = data_str_outputgap$output_gap.l1[1:obs-1], x7 = data_str_outputgap$real_interest_gap.l2[1:obs-1], x8 = data_str_outputgap$us_output_gap[1:obs-1], x9 = data_str_outputgap$rer_gap[1:obs-1], x10 = data_str_outputgap$brent_gap[1:obs-1],
                                      s_var = state_variable[2:obs], s_var_l = state_variable[1:obs-1]) 
                          ,
                          start = list(rho = par_str_ini[1], b0 = par_str_ini[2], b1 = par_str_ini[3], b2 = par_str_ini[4], b3= par_str_ini[5], b4= par_str_ini[6], b5= par_str_ini[7], b6 = par_str_ini[8], 
                                       b7 = par_str_ini[9], b8 = par_str_ini[10], b9 = par_str_ini[11], b10 = par_str_ini[12],b11 = par_str_ini[13], gamma = par_str_ini[14], c = par_str_ini[15]),
                          
                          lower = c(0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,gamma_grid[1],quantile(state_variable,probs = c_grid[1])),
                          upper = c(Inf, Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf, Inf,quantile(state_variable,probs = c_grid[length(c_grid)])),
                          control = nls.lm.control(maxiter = 1000, maxfev = 10000), nls.control(minFactor = 1/4096, warnOnly = TRUE), jac = NULL), silent = TRUE)
   
   if(try(summary(nls_str_LM),silent = TRUE)[2] %in% "try-error" || class(nls_str_LM) %in% "try-error"){next} else {
      convergence_message=if(nls_str_LM$convInfo$isConv == TRUE){convergence_message = "converged"} else {"no convergence"}
      #Thing we need to get from the model:
      # Name of the transition variable 
      name_state_variable = names_str_outputgap[i]
      # Parameters (linear and nonlinear part) 
      parameters_aux = nls_str_LM$m$getPars()
      #Std error of parameter estimates
      std_parameters_aux = summary(nls_str_LM)$coefficients[,2]
      #Transition Function (dated and ordered)
      f_estr1_aux = 1-exp(-(parameters_aux[length(parameters_aux)-1]/sd(state_variable[2:obs]))*((state_variable[2:obs]-parameters_aux[length(parameters_aux)])^2))
      #plot(f_estr1_aux, type = 'l', col = 'blue')
      #plot(state_variable[2:obs],f_estr1_aux)
      f_estr1_aux_ordered = sort(f_estr1_aux, decreasing = FALSE)
      #plot(f_lstr1_aux_ordered, col = 'red')
      #Predicted values: 
      predicted_str_aux = predict(nls_str_LM)
      #BIC (to compare against linear models) 
      BIC_str_aux = BIC(nls_str_LM)
      BIC_criteria = if(BIC_str_aux >= BIC_outputgap_nls){BIC_criteria = "Linear"} else {BIC_criteria = "No Linear"}
      #AIC (to compare against alternative STR especifications)
      AIC_str_aux = AIC(nls_str_LM)
      AIC_Criteria = if(AIC_str_aux >= AIC_outputgap_nls){AIC_Criteria = "Linear"} else {AIC_Criteria = "No Linear"}
      # Residuals of the model (with 0 at the beginning): 
      residuals_str_aux = resid(nls_str_LM)
      # SSR (necessary for evaluation stage)
      ssr_str_aux = sum((predict(nls_str_LM) - data_str_outputgap$output_gap[2:obs])^2)
      SSR_criteria = if(ssr_str_aux >= ssr_1_outputgap_nls){SSR_criteria = "Linear"} else {SSR_criteria = "No Linear"}
      #Gradient for autocorrelation test 
      gradient_aux = nls_str_LM$m$gradient()
      #Test of restriction of nonlinear parameters = 0
      F_test_nlparam = ((ssr_1_outputgap_nls-ssr_str_aux)/(summary(outputgap_eq_nls)$df[2] -summary(nls_str_LM)$df[2]))/(ssr_str_aux/(summary(nls_str_LM)$df[2]))
      pval_f_nlparam = pf(F_test_nlparam,summary(outputgap_eq_nls)$df[2]-summary(nls_str_LM)$df[2],summary(nls_str_LM)$df[2],lower.tail = FALSE)
      #Test of autocorrelation of residuals 
      #1. Orthogonalize errors: 
      ortho_error_autocorr = lm(residuals_str_aux ~ 1 +gradient_aux)
      reisduals_ortho_str_aux = resid(ortho_error_autocorr)
      #2. Matrix of lagged residuales (vt):
      vt_l_residuals = cbind(reisduals_ortho_str_aux,dplyr::lag(reisduals_ortho_str_aux,1),dplyr::lag(reisduals_ortho_str_aux,2),dplyr::lag(reisduals_ortho_str_aux,3),
                             dplyr::lag(reisduals_ortho_str_aux,4),dplyr::lag(reisduals_ortho_str_aux,5),dplyr::lag(reisduals_ortho_str_aux,6),dplyr::lag(reisduals_ortho_str_aux,7),
                             dplyr::lag(reisduals_ortho_str_aux,9),dplyr::lag(reisduals_ortho_str_aux,10),dplyr::lag(reisduals_ortho_str_aux,11),dplyr::lag(reisduals_ortho_str_aux,12))
      vt_l_residuals[is.na(vt_l_residuals)] <- 0
      vt_l_residuals = vt_l_residuals[,c(-1)]
      #3. Regression result and LM tests (F and Chi) 
      auto_corr_test = lm(residuals_str_aux ~ gradient_aux + vt_l_residuals)
      ssr_alt = sum((predict(auto_corr_test)-residuals_str_aux)^2) 
      F_corr_str = ((ssr_str_aux - ssr_alt)/(summary(nls_str_LM)$df[2]-auto_corr_test$df.residual))/(ssr_alt/(auto_corr_test$df.residual))
      pval_f_corr = pf(F_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual,auto_corr_test$df.residual,lower.tail = FALSE)
      chi_corr_str = length(residuals_str_aux)*((ssr_str_aux - ssr_alt)/ssr_str_aux) 
      pval_chi_corr = pchisq(chi_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual, lower.tail = FALSE) 
      #Result of Jarque-Bera and Shapiro Test 
      jb_test_aux = jarque.bera.test(residuals_str_aux)$p.value
      shapiro_aux = shapiro.test(residuals_str_aux)$p.value
      #Ramsey Test
      ols_str_ramsey = lm(residuals_str_aux ~ gradient_aux + predicted_str_aux + predicted_str_aux^2 + predicted_str_aux^3)
      ssr_ols_str_ramsey = sum((predict(ols_str_ramsey) - residuals_str_aux)^2)
      #F test statistic:
      f_ramsey_ols = ((ssr_str_aux - ssr_ols_str_ramsey)/(summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual))/(ssr_ols_str_ramsey/(ols_str_ramsey$df.residual))
      pval_f_str_ramsey = pf(f_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,ols_str_ramsey$df.residual, lower.tail = FALSE)
      #Chi test statistic 
      chi_ramsey_ols = length(residuals_str_aux)*((ssr_str_aux - ssr_ols_str_ramsey)/ssr_str_aux)
      pval_chi_str_ramsey = pchisq(chi_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,lower.tail = FALSE)
      #No additive Linearity Test - LM: 
      no_ad_nolin_results_direct = NULL
      for(k in eval(ncol(data_str_outputgap_model)+1):eval(ncol(data_str_outputgap_model)+length(selected_variables_direct))){
         state_2 =  data_str_outputgap[2:nrow(data_str_outputgap),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux + eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_direct = rbind(no_ad_nolin_results_direct,no_ad_nolin_results_aux)
      }
      
      no_ad_nolin_results_state = NULL
      for(k in eval(ncol(data_str_outputgap_model)+length(selected_variables_direct)+1):eval(ncol(data_str_outputgap_model)+ncol(data_str_outputgap_state))){
         state_2 =  data_str_outputgap[2:nrow(data_str_outputgap),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux +state_2+eval(state_2^2)+eval(state_2^3)+ eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_state = rbind(no_ad_nolin_results_state,no_ad_nolin_results_aux)
      }
      no_ad_nolin_results_direct = na.omit(no_ad_nolin_results_direct)
      no_ad_nolin_results_state = na.omit(no_ad_nolin_results_state)
      if(any(no_ad_nolin_results_direct[,2] <= 0.05) == TRUE){no_ad_d_nl_5_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_direct[,2] <= 0.01) == TRUE){no_ad_d_nl_1_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_1_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.05) == TRUE){no_ad_s_nl_5_pc = length(which(no_ad_nolin_results_state[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.01) == TRUE){no_ad_s_nl_1_pc = length(which(no_ad_nolin_results_state[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_1_pc="No Additive Nonlinearity"}
      #Parameter Constancy Test - LM: cannot be performed do to small sample size problem (low DF)
      #Arch test: 
      noarch_test =  ArchTest(residuals_str_aux, lags = 12)
      arch_chi_statistic = noarch_test$statistic
      arch_chi_pval = noarch_test$p.value
      
      #Final Output
      model_info = rbind(convergence_message,BIC_str_aux,BIC_criteria, AIC_str_aux, AIC_Criteria, ssr_str_aux,SSR_criteria,F_corr_str,pval_f_corr,chi_corr_str,pval_chi_corr,
                         jb_test_aux,shapiro_aux,f_ramsey_ols,pval_f_str_ramsey,chi_ramsey_ols,pval_chi_str_ramsey,no_ad_d_nl_5_pc,no_ad_d_nl_1_pc,no_ad_s_nl_5_pc,no_ad_s_nl_1_pc,arch_chi_statistic,arch_chi_pval,
                         F_test_nlparam, pval_f_nlparam) 
      colnames(model_info) = name_state_variable
      rownames(model_info) = c("Convergence", "BIC", "BIC Criteria", "AIC", "AIC Criteria", "SSR", "SSR Criteria", "F test corr", "Pval F test", "Chi2 test corr", "Pval Chi2 test",
                               "JB test", "Shapiro test", "F test Ramsey", "Pval F Ramsey", "Chi2 test Ramsey", "Pval Chi2 Ramsey", "No Ad NL D 5%","No Ad D NL 1%","No Ad NL S 5%","No Ad S NL 1%", "ARCH test stat", "ARCH test pval", 
                               "F param nl", "Pval F param nl")
      Transition_function = cbind(f_estr1_aux, f_estr1_aux_ordered, state_variable[2:obs])
      colnames(Transition_function) = c("LSTR1", "LSTR1 sorted", "State_Variable")
      final_output = list(Model_Info = model_info,Model = nls_str_LM, Transition_function = Transition_function)
      ESTR_OUTPUTGAP[[eval(i-ncol(data_str_outputgap_model))]] <- final_output
   }
}
names(ESTR_OUTPUTGAP) <- names_str_outputgap[eval(ncol(data_str_outputgap_model)+1):length(names_str_outputgap)]
proc.time() - ptm 

#Results: of model info 
ESTR_OUTPUTGAP_model_results = NULL
for(j in 1:length(ESTR_OUTPUTGAP)){aux_info = ESTR_OUTPUTGAP[[j]]$Model_Info
ESTR_OUTPUTGAP_model_results = cbind(ESTR_OUTPUTGAP_model_results,aux_info)
}
#Selected models According to BIC/AIC criteria,  
ESTR_OUTPUTGAP_selected_models = ESTR_OUTPUTGAP_model_results[,
                                                                which(ESTR_OUTPUTGAP_model_results['Convergence',] == "converged" & 
                                                                         (ESTR_OUTPUTGAP_model_results['BIC Criteria',] == "No Linear" |  ESTR_OUTPUTGAP_model_results['AIC Criteria',] == "No Linear")
                                                                      & (ESTR_OUTPUTGAP_model_results['JB test',] >= 0.05 | ESTR_OUTPUTGAP_model_results['Shapiro test',] >= 0.05)
                                                                      &  ESTR_OUTPUTGAP_model_results['Pval F param nl',] <= 0.05 &  ESTR_OUTPUTGAP_model_results['Pval F test',] >= 0.05 &  ESTR_OUTPUTGAP_model_results['Pval Chi2 test',] >= 0.05
                                                                      &  ESTR_OUTPUTGAP_model_results['Pval F Ramsey',] >= 0.05 & ESTR_OUTPUTGAP_model_results['Pval Chi2 Ramsey',] >= 0.05
                                                                      & ESTR_OUTPUTGAP_model_results['ARCH test pval',] >= 0.05 & (ESTR_OUTPUTGAP_model_results['No Ad D NL 1%',] == "No Additive Nonlinearity" | ESTR_OUTPUTGAP_model_results['No Ad D NL 1%',] < 0.25) 
                                                                      & (ESTR_OUTPUTGAP_model_results['No Ad S NL 1%',] == "No Additive Nonlinearity" | ESTR_OUTPUTGAP_model_results['No Ad S NL 1%',] < 0.25))]



LSTR2_OUTPUTGAP = list()
ptm <- proc.time()
#For Loop for estimates of LSTR1 function 
for(i in eval(ncol(data_str_outputgap_model)+1):eval(ncol(data_str_outputgap_model)+ncol(data_str_outputgap_state))){
   #Declare state variable
   state_variable = data_str_outputgap[,i]  
   #Step 1: Grid for Gamma and C parameters
   gamma_grid = seq(0.1*sd(state_variable[2:obs])^2, 32*sd(state_variable[2:obs])^2, length.out = 32) #Reescaled by transition variable 
   c1_grid = seq(0.1,0.5,0.03) #c1 parameter, from 10 th percentile to 50th percentile 
   c2_grid = seq(0.55, 0.95, 0.03) #c2 parameter, from 50th percentile to 90th percentile 
   combinations = as.data.frame(crossing(gamma = 1:length(gamma_grid), c1 = 1:length(c1_grid), c2 = 1:length(c2_grid)))
   dim(combinations)
   
   #LSTR 2 estimation
   ssr_linear_str = NULL
   parameters_linear_str = NULL
   get_init_param_o <- function(x){
      #Combinations of parameters
      gamma_cond <- gamma_grid[combinations[x,1]]
      c1_cond <-quantile(state_variable, probs = c_grid[combinations[x,2]]) #State variable transition point
      c2_cond <-quantile(state_variable, probs = c_grid[combinations[x,3]])
      # Compute transition function for the state variable
      f_lstr2 = 1/(1+exp(-(gamma_cond/sd(state_variable[2:obs])^2)*(state_variable-c1_cond)*(state_variable-c2_cond))) #This becomes a variable 
      # Estimate model by OLS (without rho term) to get initial parameters for NLS estimation
      ols_aux = lm(output_gap  ~ 1 + output_gap.l1 + real_interest_gap.l2 + rer_gap + us_output_gap + brent_gap             
                   + f_lstr2 + eval(output_gap.l1*f_lstr2) + eval(real_interest_gap.l2*f_lstr2) + eval(rer_gap*f_lstr2) + eval(us_output_gap*f_lstr2) + eval(brent_gap*f_lstr2),
                   data = cbind(data_str_outputgap,f_lstr2))
      par_ols = coef(ols_aux)
      par_ols[is.na(par_ols)]<- 0
      # Add variable to the estimation of linear model (by NLS):
      #From x11 to x15 are the STR interactions of model variables with state variables
      #From x16 to x20 are lagged interactions with rho term
      
      linear_str = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4 +b5*x5
                             -rho*b1*x6-rho*b2*x7-rho*b3*x8-rho*b4*x9-rho*b5*x10
                             +b6*x11 + b7*x12 + b8*x13 + b9*x14 + b10*x15 + b11*x16
                             - rho*b6*x17 - rho*b7*x18 - rho*b8*x19- rho*b9*x20 - rho*b10*x21 - rho*b11*x22
                             , 
                             start = list(rho = 0.03, b0 = -0.0015, b1 = 0.88, b2 = -0.023, b3= -0.084, b4= 0.23, b5= 0.004, 
                                          b6 = par_ols[7], b7 = par_ols[8], b8 = par_ols[9], b9 = par_ols[10], b10 = par_ols[11], b11 = par_ols[12]),
                             data = list(y = data_str_outputgap$output_gap[2:obs], y_l = data_str_outputgap$output_gap[1:obs-1], x1 = data_str_outputgap$output_gap.l1[2:obs],
                                         x2 = data_str_outputgap$real_interest_gap.l2[2:obs], x3 = data_str_outputgap$us_output_gap[2:obs], x4 = data_str_outputgap$rer_gap[2:obs], x5 = data_str_outputgap$brent_gap[2:obs], 
                                         x6 = data_str_outputgap$output_gap.l1[1:obs-1], x7 = data_str_outputgap$real_interest_gap.l2[1:obs-1], x8 = data_str_outputgap$us_output_gap[1:obs-1], x9 = data_str_outputgap$rer_gap[1:obs-1], x10 = data_str_outputgap$brent_gap[1:obs-1],
                                         x11 = f_lstr2[2:obs], x12 = data_str_outputgap$output_gap.l1[2:obs]*f_lstr2[2:obs], x13 = data_str_outputgap$real_interest_gap.l2[2:obs]*f_lstr2[2:obs], x14 = data_str_outputgap$us_output_gap[2:obs]*f_lstr2[2:obs], x15 = data_str_outputgap$rer_gap[2:obs]*f_lstr2[2:obs], x16 = data_str_outputgap$brent_gap[2:obs]*f_lstr2[2:obs], 
                                         x17 = f_lstr2[1:obs-1], x18 = data_str_outputgap$output_gap.l1[1:obs-1]*f_lstr2[1:obs-1], x19 = data_str_outputgap$real_interest_gap.l2[1:obs-1]*f_lstr2[1:obs-1], x20 = data_str_outputgap$us_output_gap[1:obs-1]*f_lstr2[1:obs-1], x21 = data_str_outputgap$rer_gap[1:obs-1]*f_lstr2[1:obs-1], x22 = data_str_outputgap$brent_gap[1:obs-1]*f_lstr2[1:obs-1]
                                         
                             ), 
                             trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200) , nls.control(minFactor = 1/4096, warnOnly = TRUE)), silent = FALSE)
      if(try(summary(linear_str),silent = TRUE)[2] %in% "try-error" || class(linear_str) %in% "try-error"){return(c(rep(NA, 14)))} else {
         # Get the Residual Sum of Squares of the Regression
         ssr_aux = sum((predict(linear_str) - data_str_outputgap$output_gap[2:obs]) ^ 2)  ## residual sum of squares
         # Get the model parameters
         parameters_aux = linear_str$m$getPars()
         parameters_ssr_linear_str = rbind(c(), c(parameters_aux, ssr_aux))
         return(parameters_ssr_linear_str)
      }
   }
   system.time({param_ssr_matrix = foreach(i = 1:nrow(combinations), .combine = rbind, .packages = c('doParallel', 'minpack.lm')) %dopar% {get_init_param_o(i)}})
   param_ssr_matrix = na.omit(param_ssr_matrix)
   index_selection = which(param_ssr_matrix[,ncol(param_ssr_matrix)] == min(param_ssr_matrix[,ncol(param_ssr_matrix)])) 
   par_str_ini = c(param_ssr_matrix[index_selection,-ncol(param_ssr_matrix)],gamma_grid[combinations[index_selection,1]],quantile(state_variable,probs = c1_grid[combinations[index_selection,2]]),quantile(state_variable,probs = c1_grid[combinations[index_selection,3]]))
   
   # Estimate by NLS LM
   nls_str_LM = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+b4*x4+b5*x5
                          -rho*b1*x6-rho*b2*x7-rho*b3*x8-rho*b4*x9-rho*b5*x10
                          +(1/(1+exp(-(gamma/sd(state_variable[2:obs])^2)*(s_var-c1)*(s_var-c2))))*(b6 +b7*x1+b8*x2+b9*x3+b10*x4+b11*x5)
                          -rho*(1/(1+exp(-(gamma/sd(state_variable[2:obs])^2)*(s_var_l-c1)*(s_var-c2))))*(b6+ b7*x6+b8*x7+b9*x8+b10*x9+b11*x10)
                          , 
                          data = list(y = data_str_outputgap$output_gap[2:obs], y_l = data_str_outputgap$output_gap[1:obs-1], x1 = data_str_outputgap$output_gap.l1[2:obs],
                                      x2 = data_str_outputgap$real_interest_gap.l2[2:obs], x3 = data_str_outputgap$us_output_gap[2:obs], x4 = data_str_outputgap$rer_gap[2:obs], x5 = data_str_outputgap$brent_gap[2:obs], 
                                      x6 = data_str_outputgap$output_gap.l1[1:obs-1], x7 = data_str_outputgap$real_interest_gap.l2[1:obs-1], x8 = data_str_outputgap$us_output_gap[1:obs-1], x9 = data_str_outputgap$rer_gap[1:obs-1], x10 = data_str_outputgap$brent_gap[1:obs-1],
                                      s_var = state_variable[2:obs], s_var_l = state_variable[1:obs-1]) 
                          ,
                          start = list(rho = par_str_ini[1], b0 = par_str_ini[2], b1 = par_str_ini[3], b2 = par_str_ini[4], b3= par_str_ini[5], b4= par_str_ini[6], b5= par_str_ini[7], b6 = par_str_ini[8], 
                                       b7 = par_str_ini[9], b8 = par_str_ini[10], b9 = par_str_ini[11], b10 = par_str_ini[12],b11 = par_str_ini[13], gamma = par_str_ini[14], c1 = par_str_ini[15], c2 = par_str_ini[16]),
                          
                          lower = c(0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,gamma_grid[1],quantile(state_variable,probs = c1_grid[1]),quantile(state_variable,probs = c2_grid[1])),
                          upper = c(Inf, Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf, Inf,Inf,Inf,quantile(state_variable,probs = c1_grid[length(c1_grid)]),quantile(state_variable,probs = c2_grid[length(c2_grid)])),
                          control = nls.lm.control(maxiter = 200), nls.control(minFactor = 1/4096, warnOnly = TRUE), jac = NULL), silent = TRUE)
   if(try(summary(nls_str_LM),silent = TRUE)[2] %in% "try-error" || class(nls_str_LM) %in% "try-error"){next} else {
      convergence_message=if(nls_str_LM$convInfo$isConv == TRUE){convergence_message = "converged"} else {"no convergence"}
      #Thing we need to get from the model:
      # Name of the transition variable 
      name_state_variable = names_str_outputgap[i]
      # Parameters (linear and nonlinear part) 
      parameters_aux = nls_str_LM$m$getPars()
      #Std error of parameter estimates
      std_parameters_aux = summary(nls_str_LM)$coefficients[,2]
      #Transition Function (dated and ordered)
      f_lstr2_aux = 1/(1+exp(-(parameters_aux[length(parameters_aux)-2]/sd(state_variable[2:obs])^2)*(state_variable[2:obs]-parameters_aux[length(parameters_aux)-1])*(state_variable[2:obs]-parameters_aux[length(parameters_aux)])))
      #plot(f_lstr2_aux, type = 'l', col = 'blue')
      #plot(state_variable[2:obs],f_lstr2_aux)
      f_lstr2_aux_ordered = sort(f_lstr2_aux, decreasing = FALSE)
      #plot(f_lstr2_aux_ordered, col = 'red')
      #Predicted values: 
      predicted_str_aux = predict(nls_str_LM)
      #BIC (to compare against linear models) 
      BIC_str_aux = BIC(nls_str_LM)
      BIC_criteria = if(BIC_str_aux >= BIC_outputgap_nls){BIC_criteria = "Linear"} else {BIC_criteria = "No Linear"}
      #AIC (to compare against alternative STR especifications)
      AIC_str_aux = AIC(nls_str_LM)
      AIC_Criteria = if(AIC_str_aux >= AIC_outputgap_nls){AIC_Criteria = "Linear"} else {AIC_Criteria = "No Linear"}
      # Residuals of the model (with 0 at the beginning): 
      residuals_str_aux = resid(nls_str_LM)
      # SSR (necessary for evaluation stage)
      ssr_str_aux = sum((predict(nls_str_LM) - data_str_outputgap$output_gap[2:obs])^2)
      SSR_criteria = if(ssr_str_aux >= ssr_1_outputgap_nls){SSR_criteria = "Linear"} else {SSR_criteria = "No Linear"}
      #Gradient for autocorrelation test 
      gradient_aux = nls_str_LM$m$gradient()
      #Test of restriction of nonlinear parameters = 0
      F_test_nlparam = ((ssr_1_outputgap_nls-ssr_str_aux)/(summary(outputgap_eq_nls)$df[2] -summary(nls_str_LM)$df[2]))/(ssr_str_aux/(summary(nls_str_LM)$df[2]))
      pval_f_nlparam = pf(F_test_nlparam,summary(outputgap_eq_nls)$df[2]-summary(nls_str_LM)$df[2],summary(nls_str_LM)$df[2],lower.tail = FALSE)
      #Test of autocorrelation of residuals 
      #1. Orthogonalize errors: 
      ortho_error_autocorr = lm(residuals_str_aux ~ 1 +gradient_aux)
      reisduals_ortho_str_aux = resid(ortho_error_autocorr)
      #2. Matrix of lagged residuales (vt):
      vt_l_residuals = cbind(reisduals_ortho_str_aux,dplyr::lag(reisduals_ortho_str_aux,1),dplyr::lag(reisduals_ortho_str_aux,2),dplyr::lag(reisduals_ortho_str_aux,3),
                             dplyr::lag(reisduals_ortho_str_aux,4),dplyr::lag(reisduals_ortho_str_aux,5),dplyr::lag(reisduals_ortho_str_aux,6),dplyr::lag(reisduals_ortho_str_aux,7),
                             dplyr::lag(reisduals_ortho_str_aux,9),dplyr::lag(reisduals_ortho_str_aux,10),dplyr::lag(reisduals_ortho_str_aux,11),dplyr::lag(reisduals_ortho_str_aux,12))
      vt_l_residuals[is.na(vt_l_residuals)] <- 0
      vt_l_residuals = vt_l_residuals[,c(-1)]
      #3. Regression result and LM tests (F and Chi) 
      auto_corr_test = lm(residuals_str_aux ~ gradient_aux + vt_l_residuals)
      ssr_alt = sum((predict(auto_corr_test)-residuals_str_aux)^2) 
      F_corr_str = ((ssr_str_aux - ssr_alt)/(summary(nls_str_LM)$df[2]-auto_corr_test$df.residual))/(ssr_alt/(auto_corr_test$df.residual))
      pval_f_corr = pf(F_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual,auto_corr_test$df.residual,lower.tail = FALSE)
      chi_corr_str = length(residuals_str_aux)*((ssr_str_aux - ssr_alt)/ssr_str_aux) 
      pval_chi_corr = pchisq(chi_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual, lower.tail = FALSE) 
      #Result of Jarque-Bera and Shapiro Test 
      jb_test_aux = jarque.bera.test(residuals_str_aux)$p.value
      shapiro_aux = shapiro.test(residuals_str_aux)$p.value
      #Ramsey Test
      ols_str_ramsey = lm(residuals_str_aux ~ gradient_aux + predicted_str_aux + predicted_str_aux^2 + predicted_str_aux^3)
      ssr_ols_str_ramsey = sum((predict(ols_str_ramsey) - residuals_str_aux)^2)
      #F test statistic:
      f_ramsey_ols = ((ssr_str_aux - ssr_ols_str_ramsey)/(summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual))/(ssr_ols_str_ramsey/(ols_str_ramsey$df.residual))
      pval_f_str_ramsey = pf(f_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,ols_str_ramsey$df.residual, lower.tail = FALSE)
      #Chi test statistic 
      chi_ramsey_ols = length(residuals_str_aux)*((ssr_str_aux - ssr_ols_str_ramsey)/ssr_str_aux)
      pval_chi_str_ramsey = pchisq(chi_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,lower.tail = FALSE)
      #No additive Linearity Test - LM: 
      no_ad_nolin_results_direct = NULL
      for(k in eval(ncol(data_str_outputgap_model)+1):eval(ncol(data_str_outputgap_model)+length(selected_variables_direct))){
         state_2 =  data_str_outputgap[2:nrow(data_str_outputgap),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux + eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_direct = rbind(no_ad_nolin_results_direct,no_ad_nolin_results_aux)
      }
      
      no_ad_nolin_results_state = NULL
      for(k in eval(ncol(data_str_outputgap_model)+length(selected_variables_direct)+1):eval(ncol(data_str_outputgap_model)+ncol(data_str_outputgap_state))){
         state_2 =  data_str_outputgap[2:nrow(data_str_outputgap),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux +state_2+eval(state_2^2)+eval(state_2^3)+ eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_state = rbind(no_ad_nolin_results_state,no_ad_nolin_results_aux)
      }
      no_ad_nolin_results_direct = na.omit(no_ad_nolin_results_direct)
      no_ad_nolin_results_state = na.omit(no_ad_nolin_results_state)
      if(any(no_ad_nolin_results_direct[,2] <= 0.05) == TRUE){no_ad_d_nl_5_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_direct[,2] <= 0.01) == TRUE){no_ad_d_nl_1_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_1_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.05) == TRUE){no_ad_s_nl_5_pc = length(which(no_ad_nolin_results_state[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.01) == TRUE){no_ad_s_nl_1_pc = length(which(no_ad_nolin_results_state[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_1_pc="No Additive Nonlinearity"}
      #Parameter Constancy Test - LM: cannot be performed do to small sample size problem (low DF)
      #Arch test: 
      noarch_test =  ArchTest(residuals_str_aux, lags = 12)
      arch_chi_statistic = noarch_test$statistic
      arch_chi_pval = noarch_test$p.value
      
      #Final Output
      model_info = rbind(convergence_message,BIC_str_aux,BIC_criteria, AIC_str_aux, AIC_Criteria, ssr_str_aux,SSR_criteria,F_corr_str,pval_f_corr,chi_corr_str,pval_chi_corr,
                         jb_test_aux,shapiro_aux,f_ramsey_ols,pval_f_str_ramsey,chi_ramsey_ols,pval_chi_str_ramsey,no_ad_d_nl_5_pc,no_ad_d_nl_1_pc,no_ad_s_nl_5_pc,no_ad_s_nl_1_pc,arch_chi_statistic,arch_chi_pval, 
                         F_test_nlparam, pval_f_nlparam) 
      colnames(model_info) = name_state_variable
      rownames(model_info) = c("Convergence", "BIC", "BIC Criteria", "AIC", "AIC Criteria", "SSR", "SSR Criteria", "F test corr", "Pval F test", "Chi2 test corr", "Pval Chi2 test",
                               "JB test", "Shapiro test", "F test Ramsey", "Pval F Ramsey", "Chi2 test Ramsey", "Pval Chi2 Ramsey", "No Ad NL D 5%","No Ad D NL 1%","No Ad NL S 5%","No Ad S NL 1%", "ARCH test stat", "ARCH test pval", 
                               "F test nl", "Pval F nl")
      Transition_function = cbind(f_lstr2_aux, f_lstr2_aux_ordered, state_variable[2:obs])
      colnames(Transition_function) = c("LSTR2", "LSTR2 sorted", "State_Variable")
      final_output = list(Model_Info = model_info,Model = nls_str_LM, Transition_function = Transition_function)
      LSTR2_OUTPUTGAP[[eval(i-ncol(data_str_outputgap_model))]] <- final_output
   }
}
names(LSTR2_OUTPUTGAP) <- names_str_inf[eval(ncol(data_str_outputgap_model)+1):length(names_str_outputgap)]
proc.time() - ptm 

#Results: of model info 
LSTR2_OUTPUTGAP_model_results = NULL
for(j in 1:length(LSTR2_OUTPUTGAP)){aux_info = LSTR2_OUTPUTGAP[[j]]$Model_Info
LSTR2_OUTPUTGAP_model_results = cbind(LSTR2_OUTPUTGAP_model_results,aux_info)
}
#Selected models According to BIC/AIC criteria,  
LSTR2_OUTPUTGAP_selected_models = LSTR2_OUTPUTGAP_model_results[,
                                                              which(LSTR2_OUTPUTGAP_model_results['Convergence',] == "converged" & 
                                                                       (LSTR2_OUTPUTGAP_model_results['BIC Criteria',] == "No Linear" |  LSTR2_OUTPUTGAP_model_results['AIC Criteria',] == "No Linear")
                                                                    & (LSTR2_OUTPUTGAP_model_results['JB test',] >= 0.05 | LSTR2_OUTPUTGAP_model_results['Shapiro test',] >= 0.05)
                                                                    &  LSTR2_OUTPUTGAP_model_results['Pval F nl',] <= 0.05 &  LSTR2_OUTPUTGAP_model_results['Pval F test',] >= 0.05 &  LSTR2_OUTPUTGAP_model_results['Pval Chi2 test',] >= 0.05
                                                                    &  LSTR2_OUTPUTGAP_model_results['Pval F Ramsey',] >= 0.05 & LSTR2_OUTPUTGAP_model_results['Pval Chi2 Ramsey',] >= 0.05
                                                                    & LSTR2_OUTPUTGAP_model_results['ARCH test pval',] >= 0.05 & (LSTR2_OUTPUTGAP_model_results['No Ad D NL 1%',] == "No Additive Nonlinearity" | LSTR2_OUTPUTGAP_model_results['No Ad D NL 1%',] < 0.25) 
                                                                    & (LSTR2_OUTPUTGAP_model_results['No Ad S NL 1%',] == "No Additive Nonlinearity" | LSTR2_OUTPUTGAP_model_results['No Ad S NL 1%',] < 0.25))]

# ------- Interest Rate Equation ------------- #
# Selection of variables where LM test detected non-linearity
names_interest_direct = colnames(LM_interest_direct)
names_interest_state = colnames(LM_interest_state)
selected_variables_direct = names_interest_direct[which(LM_interest_direct[c("LM result"),] == "Nonlinear")]
selected_variables_state = names_interest_state[which(LM_interest_state[c("LM result"),] == "Nonlinear")]
selected_variables = c(selected_variables_direct,selected_variables_state)

data_str_interest_state = data_state_interest[,selected_variables] 
data_str_interest_state = xts(data_str_interest_state, index_dates_data_state_final)
data_str_interest_model = as.data.frame(cbind(data_state_final$interest_rate, data_state_final$resid_interest, data_state_final$inflation_gap, data_state_final$output_gap,data_state_final$interest_rate.l1))
data_str_interest_model = xts(data_str_interest_model, index_dates_data_state_final)
data_str_interest = cbind.xts(data_str_interest_model, data_str_interest_state)

obs = nrow(data_str_interest)
data_str_interest = as.data.frame(data_str_interest)
names_str_interest = colnames(data_str_interest)

#Steps (for each nonlinear or state variable)
#1. Calculate Grid for Gamma and C parameters
#2. USe the Grid to compute the transition function (logistic/exponential)
#3. Add transition function to the model and estimate model by OLS - NLS
#4. Obtain SSR for each iteration, and choose the parameters that minimize it
#5. Use parameters in 4 as starting values to NLS estimate of STR

LSTR1_INTEREST = list()
ptm <- proc.time()
#For Loop for estimates of LSTR1 function 
for(i in eval(ncol(data_str_interest_model)+1):eval(ncol(data_str_interest_model)+ncol(data_str_interest_state))){
   #Declare state variable
   state_variable = data_str_interest[,i]  
   #Grid for Gamma and C parameters
   gamma_grid = seq(0.1*sd(state_variable[2:obs]), 32*sd(state_variable[2:obs]), length.out = 32) #Reescaled by transition variable 
   c_grid = seq(0.1, 0.95, 0.02) #From 10th percentile to 90th percentile
   combinations = as.data.frame(crossing(gamma = 1:length(gamma_grid), c = 1:length(c_grid)))
   dim(combinations)
   
   #LSTR 1 estimation
   ssr_linear_str = NULL
   parameters_linear_str = NULL
   get_init_param_int <- function(x){
      #Combinations of parameters
      gamma_cond <- gamma_grid[combinations[x,1]]
      c_cond <-quantile(state_variable[2:obs], probs = c_grid[combinations[x,2]]) #State variable transition point
      # Compute transition function for the state variable
      f_lstr1 = 1/(1+exp(-(gamma_cond/sd(state_variable[2:obs]))*(state_variable-c_cond))) #This becomes a variable 
      # Estimate model by OLS (without rho term) to get initial parameters for NLS estimation
      ols_aux = lm(interest_rate ~ 1 + interest_rate.l1+inflation_gap + output_gap  
                   +f_lstr1+ eval(inflation_gap*f_lstr1) + eval(output_gap*f_lstr1) + eval(interest_rate.l1*f_lstr1),
                   data = cbind(data_str_interest,f_lstr1))
      par_ols = coef(ols_aux)
      par_ols[is.na(par_ols)]<-0
      # Add variable to the estimation of linear model (by NLS):
      #From x7 to x9 are the STR interactions of model variables with state variables
      #From x11 to x12 are lagged interactions with rho term
      linear_str = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+
                                -rho*b1*x4-rho*b2*x5-rho*b3*x6
                             +b4*x7 + b5*x8 + b6*x9 + b7*x10
                             -rho*b4*x11 - rho*b5*x12 - rho*b6*x13 - rho*b7*x14
                             , 
                             start = list(rho = 0.03, b0 = 0.0, b1 = 0.9, b2 = 0.4, b3= 0.15, 
                                          b4 = par_ols[5], b5 = par_ols[6], b6 = par_ols[7], b7 = par_ols[8]),
                             data = list(y = data_str_interest$interest_rate[2:obs], y_l = data_str_interest$interest_rate[1:obs-1], x1 = data_str_interest$interest_rate.l1[2:obs],
                                         x2 = data_str_interest$inflation_gap[2:obs], x3 = data_str_interest$output_gap[2:obs], x4 = data_str_interest$interest_rate.l1[1:obs-1], 
                                         x5 = data_str_interest$inflation_gap[1:obs-1], x6 = data_str_interest$output_gap[1:obs-1],
                                         x7 = f_lstr1[2:obs], x8 = data_str_interest$interest_rate.l1[2:obs]*f_lstr1[2:obs],x9 = data_str_interest$inflation_gap[2:obs]*f_lstr1[2:obs], x10 = data_str_interest$output_gap[2:obs]*f_lstr1[2:obs],
                                         x11 = f_lstr1[1:obs-1], x12= data_str_interest$interest_rate.l1[1:obs-1]*f_lstr1[1:obs-1],x13 = data_str_interest$inflation_gap[1:obs-1]*f_lstr1[1:obs-1], x14 = data_str_interest$output_gap[1:obs-1]*f_lstr1[1:obs-1]),
                             trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200), nls.control(minFactor = 1/4096, warnOnly = TRUE)), silent = TRUE)
      if(try(summary(linear_str), silent = TRUE)[2] %in% "try-error"||class(linear_str) %in% "try-error"){return(c(rep(NA,10)))} else {
         # Get the Residual Sum of Squares of the Regression
         ssr_aux = sum((predict(linear_str) - data_str_interest$interest_rate[2:obs]) ^ 2)  ## residual sum of squares
         # Get the model parameters
         parameters_aux = linear_str$m$getPars()
         parameters_ssr_linear_str = rbind(c(), c(parameters_aux, ssr_aux))
         return(parameters_ssr_linear_str)
      }
   }
   system.time({param_ssr_matrix = foreach(i = 1:nrow(combinations), .combine = rbind, .packages = c('doParallel', 'minpack.lm')) %dopar% {get_init_param_int(i)}})
   param_ssr_matrix = na.omit(param_ssr_matrix)
   index_selection = which(param_ssr_matrix[,ncol(param_ssr_matrix)] == min(param_ssr_matrix[,ncol(param_ssr_matrix)])) 
   par_str_ini = c(param_ssr_matrix[index_selection,-ncol(param_ssr_matrix)],gamma_grid[combinations[index_selection,1]],quantile(state_variable,probs = c_grid[combinations[index_selection,2]]))

   # Estimate by NLS LM
   nls_str_LM = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+
                             -rho*b1*x4-rho*b2*x5-rho*b3*x6
                          +(1/(1+exp(-(gamma/sd(state_variable[2:obs]))*(s_var-c))))*(b4 + b5*x1 + b6*x2 + b7*x3)
                          -rho*(1/(1+exp(-(gamma/sd(state_variable[2:obs]))*(s_var_l-c))))*(b4 + b5*x4 + b6*x5 + b7*x6)
                          , 
                          data = list(y = data_str_interest$interest_rate[2:obs], y_l = data_str_interest$interest_rate[1:obs-1], x1 = data_str_interest$interest_rate.l1[2:obs],
                                      x2 = data_str_interest$inflation_gap[2:obs], x3 = data_str_interest$output_gap[2:obs], x4 = data_str_interest$interest_rate.l1[1:obs-1], 
                                      x5 = data_str_interest$inflation_gap[1:obs-1], x6 = data_str_interest$output_gap[1:obs-1],
                          s_var = state_variable[2:obs], s_var_l = state_variable[1:obs-1]),
                          start = list(rho = par_str_ini[1], b0 = par_str_ini[2], b1 = par_str_ini[3], b2 = par_str_ini[4], b3= par_str_ini[5], b4= par_str_ini[6], b5= par_str_ini[7], b6 = par_str_ini[8], b7 = par_str_ini[9], 
                                       gamma = par_str_ini[10], c = par_str_ini[11]),
                          
                          lower = c(0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,gamma_grid[1],quantile(state_variable,probs = c_grid[1])),
                          upper = c(Inf, Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,quantile(state_variable,probs = c_grid[length(c_grid)])),
                          control = nls.lm.control(maxiter = 200),nls.control(minFactor = 1/4096, warnOnly = TRUE), jac = NULL), silent = TRUE)
   if(try(summary(nls_str_LM),silent = TRUE)[2] %in% "try-error" || class(nls_str_LM) %in% "try-error"){next} else {
      convergence_message=if(nls_str_LM$convInfo$isConv == TRUE){convergence_message = "converged"} else {"no convergence"}
      #Thing we need to get from the model:
      # Name of the transition variable 
      name_state_variable = names_str_interest[i]
      # Parameters (linear and nonlinear part) 
      parameters_aux = nls_str_LM$m$getPars()
      #Std error of parameter estimates
      std_parameters_aux = summary(nls_str_LM)$coefficients[,2]
      #Transition Function (dated and ordered)
      f_lstr1_aux = 1/(1+exp(-(parameters_aux[length(parameters_aux)-1]/sd(state_variable[2:obs]))*(state_variable[2:obs]-parameters_aux[length(parameters_aux)])))
      #plot(state_variable[2:obs],f_lstr1_aux, col = 'blue')
      f_lstr1_aux_ordered = sort(f_lstr1_aux, decreasing = FALSE)
       #plot(f_lstr1_aux_ordered, col = 'red')
      #Predicted values: 
      predicted_str_aux = predict(nls_str_LM)
      #BIC (to compare against linear models) 
      BIC_str_aux = BIC(nls_str_LM)
      BIC_criteria = if(BIC_str_aux >= BIC_interest_nls){BIC_criteria = "Linear"} else {BIC_criteria = "No Linear"}
      #AIC (to compare against alternative STR especifications)
      AIC_str_aux = AIC(nls_str_LM)
      AIC_Criteria = if(AIC_str_aux >= AIC_interest_nls){AIC_Criteria = "Linear"} else {AIC_Criteria = "No Linear"}
      # Residuals of the model (with 0 at the beginning): 
      residuals_str_aux = resid(nls_str_LM)
      # SSR (necessary for evaluation stage)
      ssr_str_aux = sum((predict(nls_str_LM) -data_str_interest$interest_rate[2:obs])^2)
      SSR_criteria = if(ssr_str_aux >= ssr_1_interest_nls){SSR_criteria = "Linear"} else {SSR_criteria = "No Linear"}
      #Gradient for autocorrelation test 
      gradient_aux = nls_str_LM$m$gradient()
      #Test of restriction of nonlinear parameters = 0
      F_test_nlparam = ((ssr_1_interest_nls-ssr_str_aux)/(summary(interest_eq_nls)$df[2] -summary(nls_str_LM)$df[2]))/(ssr_str_aux/(summary(nls_str_LM)$df[2]))
      pval_f_nlparam = pf(F_test_nlparam,summary(interest_eq_nls)$df[2]-summary(nls_str_LM)$df[2],summary(nls_str_LM)$df[2],lower.tail = FALSE)
      #Test of autocorrelation of residuals 
      #1. Orthogonalize errors: 
      ortho_error_autocorr = lm(residuals_str_aux ~ 1+gradient_aux)
      reisduals_ortho_str_aux = resid(ortho_error_autocorr)
      #2. Matrix of lagged residuales (vt):
      vt_l_residuals = cbind(reisduals_ortho_str_aux,dplyr::lag(reisduals_ortho_str_aux,1),dplyr::lag(reisduals_ortho_str_aux,2),dplyr::lag(reisduals_ortho_str_aux,3),
                             dplyr::lag(reisduals_ortho_str_aux,4),dplyr::lag(reisduals_ortho_str_aux,5),dplyr::lag(reisduals_ortho_str_aux,6),dplyr::lag(reisduals_ortho_str_aux,7),
                             dplyr::lag(reisduals_ortho_str_aux,9),dplyr::lag(reisduals_ortho_str_aux,10),dplyr::lag(reisduals_ortho_str_aux,11),dplyr::lag(reisduals_ortho_str_aux,12))
      vt_l_residuals[is.na(vt_l_residuals)] <- 0
      vt_l_residuals = vt_l_residuals[,c(-1)]
      #3. Regression result and LM tests (F and Chi) 
      auto_corr_test = lm(residuals_str_aux ~ gradient_aux + vt_l_residuals)
      ssr_alt = sum((predict(auto_corr_test)-residuals_str_aux)^2) 
      F_corr_str = ((ssr_str_aux - ssr_alt)/(summary(nls_str_LM)$df[2]-auto_corr_test$df.residual))/(ssr_alt/(auto_corr_test$df.residual))
      pval_f_corr = pf(F_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual,auto_corr_test$df.residual,lower.tail = FALSE)
      chi_corr_str = length(residuals_str_aux)*((ssr_str_aux - ssr_alt)/ssr_str_aux) 
      pval_chi_corr = pchisq(chi_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual, lower.tail = FALSE) 
      #Result of Jarque-Bera and Shapiro Test 
      jb_test_aux = jarque.bera.test(residuals_str_aux)$p.value
      shapiro_aux = shapiro.test(residuals_str_aux)$p.value
      #Ramsey Test
      ols_str_ramsey = lm(residuals_str_aux ~ gradient_aux + predicted_str_aux + predicted_str_aux^2 + predicted_str_aux^3)
      ssr_ols_str_ramsey = sum((predict(ols_str_ramsey) - residuals_str_aux)^2)
      #F test statistic:
      f_ramsey_ols = ((ssr_str_aux - ssr_ols_str_ramsey)/(summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual))/(ssr_ols_str_ramsey/(ols_str_ramsey$df.residual))
      pval_f_str_ramsey = pf(f_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,ols_str_ramsey$df.residual, lower.tail = FALSE)
      #Chi test statistic 
      chi_ramsey_ols = length(residuals_str_aux)*(ssr_str_aux - ssr_ols_str_ramsey)/ssr_str_aux
      pval_chi_str_ramsey = pchisq(chi_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,lower.tail = FALSE)
      #No additive Linearity Test - LM: 
      no_ad_nolin_results_direct = NULL
      for(k in eval(ncol(data_str_interest_model)+1):eval(ncol(data_str_interest_model)+length(selected_variables_direct))){
         state_2 =  data_str_inf[2:nrow(data_str_interest),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux + eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_direct = rbind(no_ad_nolin_results_direct,no_ad_nolin_results_aux)
      }
      
      no_ad_nolin_results_state = NULL
      for(k in eval(ncol(data_str_interest_model)+length(selected_variables_direct)+1):eval(ncol(data_str_interest_model)+ncol(data_str_interest_state))){
         state_2 =  data_str_interest[2:nrow(data_str_interest),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux +state_2+eval(state_2^2)+eval(state_2^3)+ eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_state = rbind(no_ad_nolin_results_state,no_ad_nolin_results_aux)
      }
      no_ad_nolin_results_direct = na.omit(no_ad_nolin_results_direct)
      no_ad_nolin_results_state = na.omit(no_ad_nolin_results_state)
      if(any(no_ad_nolin_results_direct[,2] <= 0.05) == TRUE){no_ad_d_nl_5_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_direct[,2] <= 0.01) == TRUE){no_ad_d_nl_1_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_1_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.05) == TRUE){no_ad_s_nl_5_pc = length(which(no_ad_nolin_results_state[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.01) == TRUE){no_ad_s_nl_1_pc = length(which(no_ad_nolin_results_state[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_1_pc="No Additive Nonlinearity"}
      #Parameter Constancy Test - LM: cannot be performed do to small sample size problem (low DF)
      #Arch test: 
      noarch_test =  ArchTest(residuals_str_aux, lags = 12)
      arch_chi_statistic = noarch_test$statistic
      arch_chi_pval = noarch_test$p.value
      
      #Final Output
      model_info = rbind(convergence_message,BIC_str_aux,BIC_criteria, AIC_str_aux, AIC_Criteria, ssr_str_aux,SSR_criteria,F_corr_str,pval_f_corr,chi_corr_str,pval_chi_corr,
                         jb_test_aux,shapiro_aux,f_ramsey_ols,pval_f_str_ramsey,chi_ramsey_ols,pval_chi_str_ramsey,no_ad_d_nl_5_pc,no_ad_d_nl_1_pc,no_ad_s_nl_5_pc,no_ad_s_nl_1_pc,arch_chi_statistic,arch_chi_pval,
                         F_test_nlparam, pval_f_nlparam) 
      colnames(model_info) = name_state_variable
      rownames(model_info) = c("Convergence", "BIC", "BIC Criteria", "AIC", "AIC Criteria", "SSR", "SSR Criteria", "F test corr", "Pval F test", "Chi2 test corr", "Pval Chi2 test",
                               "JB test", "Shapiro test", "F test Ramsey", "Pval F Ramsey", "Chi2 test Ramsey", "Pval Chi2 Ramsey", "No Ad NL D 5%","No Ad D NL 1%","No Ad NL S 5%","No Ad S NL 1%", "ARCH test stat", "ARCH test pval", 
                               "F test nl", "Pval F nl")
      Transition_function = cbind(f_lstr1_aux, f_lstr1_aux_ordered, state_variable[2:obs])
      colnames(Transition_function) = c("LSTR1", "LSTR1 sorted", "State_Variable")
      final_output = list(Model_Info = model_info,Model = nls_str_LM, Transition_function = Transition_function)
      LSTR1_INTEREST[[eval(i-ncol(data_str_interest_model))]] <- final_output
   }
}
names(LSTR1_INTEREST) <- names_str_inf[eval(ncol(data_str_interest_model)+1):length(names_str_interest)]
proc.time() - ptm 

#Results: of model info 
LSTR1_INTEREST_model_results = NULL
for(j in 1:length(LSTR1_INTEREST)){aux_info = LSTR1_INTEREST[[j]]$Model_Info
LSTR1_INTEREST_model_results = cbind(LSTR1_INTEREST_model_results,aux_info)
}
#Selected models According to BIC/AIC criteria,  
LSTR1_INTEREST_selected_models = LSTR1_INTEREST_model_results[,
                                                                which(LSTR1_INTEREST_model_results['Convergence',] == "converged" & 
                                                                         (LSTR1_INTEREST_model_results['BIC Criteria',] == "No Linear" |  LSTR1_INTEREST_model_results['AIC Criteria',] == "No Linear")
                                                                      & (LSTR1_INTEREST_model_results['JB test',] >= 0.05 | LSTR1_INTEREST_model_results['Shapiro test',] >= 0.05)
                                                                      &  LSTR1_INTEREST_model_results['Pval F nl',] <= 0.05 &  LSTR1_INTEREST_model_results['Pval F test',] >= 0.05 &  LSTR1_INTEREST_model_results['Pval Chi2 test',] >= 0.05
                                                                      &  LSTR1_INTEREST_model_results['Pval F Ramsey',] >= 0.05 & LSTR1_INTEREST_model_results['Pval Chi2 Ramsey',] >= 0.05
                                                                      & LSTR1_INTEREST_model_results['ARCH test pval',] >= 0.05 & (LSTR1_INTEREST_model_results['No Ad D NL 1%',] == "No Additive Nonlinearity" | LSTR1_INTEREST_model_results['No Ad D NL 1%',] < 0.25) 
                                                                      & (LSTR1_INTEREST_model_results['No Ad S NL 1%',] == "No Additive Nonlinearity" | LSTR1_INTEREST_model_results['No Ad S NL 1%',] < 0.25))]

write.csv(LSTR1_INTEREST_selected_models, "LSTR1_INTEREST_selected.csv")


# Selection of variables where LM test detected non-linearity and the Wald test suggested LSTR2-ESTR model specification:
names_interest_direct = colnames(LM_interest_direct)
names_interest_state = colnames(LM_interest_state)
selected_variables_direct = names_interest_direct[which(LM_interest_direct[c("LM result"),] == "Nonlinear" & LM_interest_direct[c("STR model F"),] == "LSTR2-ESTR")]
selected_variables_state = names_interest_state[which(LM_interest_state[c("LM result"),] == "Nonlinear" & LM_interest_state[c("STR model F"),] == "LSTR2-ESTR")]
selected_variables = c(selected_variables_direct,selected_variables_state)

data_str_interest_state = data_state_interest[,selected_variables] 
data_str_interest_state = xts(data_str_interest_state, index_dates_data_state_final)
data_str_interest_model = as.data.frame(cbind(data_state_final$interest_rate, data_state_final$resid_interest, data_state_final$inflation_gap, data_state_final$output_gap,data_state_final$interest_rate.l1))
data_str_interest_model = xts(data_str_interest_model, index_dates_data_state_final)
data_str_interest = cbind.xts(data_str_interest_model, data_str_interest_state)

obs = nrow(data_str_interest)
data_str_interest = as.data.frame(data_str_interest)
names_str_interest = colnames(data_str_interest)


ESTR_INTEREST = list()
ptm <- proc.time()
#For Loop for estimates of LSTR1 function 
for(i in eval(ncol(data_str_interest_model)+1):eval(ncol(data_str_interest_model)+ncol(data_str_interest_state))){
   #Declare state variable
   state_variable = data_str_interest[,i]  
   #Grid for Gamma and C parameters
   gamma_grid = seq(0.1*sd(state_variable[2:obs]), 32*sd(state_variable[2:obs]), length.out = 32) #Reescaled by transition variable 
   c_grid = seq(0.1, 0.95, 0.02) #From 10th percentile to 90th percentile
   combinations = as.data.frame(crossing(gamma = 1:length(gamma_grid), c = 1:length(c_grid)))
   dim(combinations)
   
   #ESTR 1 estimation
   ssr_linear_str = NULL
   parameters_linear_str = NULL
   get_init_param_int <- function(x){
      #Combinations of parameters
      gamma_cond <- gamma_grid[combinations[x,1]]
      c_cond <-quantile(state_variable[2:obs], probs = c_grid[combinations[x,2]]) #State variable transition point
      # Compute transition function for the state variable
      f_estr1 = 1-exp(-(gamma_cond/sd(state_variable[2:obs]))*((state_variable-c_cond)^2)) #This becomes a variable 
      # Estimate model by OLS (without rho term) to get initial parameters for NLS estimation
      ols_aux = lm(interest_rate ~ 1 + interest_rate.l1+inflation_gap + output_gap  
                   +f_estr1+ eval(inflation_gap*f_estr1) + eval(output_gap*f_estr1) + eval(interest_rate.l1*f_estr1),
                   data = cbind(data_str_interest,f_estr1))
      par_ols = coef(ols_aux)
      par_ols[is.na(par_ols)]<- 0
      # Add variable to the estimation of linear model (by NLS):
      #From x7 to x9 are the STR interactions of model variables with state variables
      #From x11 to x12 are lagged interactions with rho term
      linear_str = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+
                                -rho*b1*x4-rho*b2*x5-rho*b3*x6
                             +b4*x7 + b5*x8 + b6*x9 + b7*x10
                             -rho*b4*x11 - rho*b5*x12 - rho*b6*x13 - rho*b7*x14
                             , 
                             start = list(rho = 0.03, b0 = 0.0, b1 = 0.9, b2 = 0.4, b3= 0.15, 
                                          b4 = par_ols[5], b5 = par_ols[6], b6 = par_ols[7], b7 = par_ols[8]),
                             data = list(y = data_str_interest$interest_rate[2:obs], y_l = data_str_interest$interest_rate[1:obs-1], x1 = data_str_interest$interest_rate.l1[2:obs],
                                         x2 = data_str_interest$inflation_gap[2:obs], x3 = data_str_interest$output_gap[2:obs], x4 = data_str_interest$interest_rate.l1[1:obs-1], 
                                         x5 = data_str_interest$inflation_gap[1:obs-1], x6 = data_str_interest$output_gap[1:obs-1],
                                         x7 = f_estr1[2:obs], x8 = data_str_interest$interest_rate.l1[2:obs]*f_estr1[2:obs],x9 = data_str_interest$inflation_gap[2:obs]*f_estr1[2:obs], x10 = data_str_interest$output_gap[2:obs]*f_estr1[2:obs],
                                         x11 =f_estr1[1:obs-1], x12 = data_str_interest$interest_rate.l1[1:obs-1]*f_estr1[1:obs-1],x13 = data_str_interest$inflation_gap[1:obs-1]*f_estr1[1:obs-1], x14 = data_str_interest$output_gap[1:obs-1]*f_estr1[1:obs-1]),
                             trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200), nls.control(minFactor = 1/4096, warnOnly = TRUE)), silent = TRUE)
      if(try(summary(linear_str), silent = TRUE)[2] %in% "try-error"||class(linear_str) %in% "try-error"){return(c(rep(NA, 10)))} else {
         # Get the Residual Sum of Squares of the Regression
         ssr_aux = sum((predict(linear_str) - data_str_interest$interest_rate[2:obs]) ^ 2)  ## residual sum of squares
         # Get the model parameters
         parameters_aux = linear_str$m$getPars()
         parameters_ssr_linear_str = rbind(c(), c(parameters_aux, ssr_aux))
         return(parameters_ssr_linear_str)
      }
   }
   system.time({param_ssr_matrix = foreach(i = 1:nrow(combinations), .combine = rbind, .packages = c('doParallel', 'minpack.lm')) %dopar% {get_init_param_int(i)}})
   param_ssr_matrix = na.omit(param_ssr_matrix)
   index_selection = which(param_ssr_matrix[,ncol(param_ssr_matrix)] == min(param_ssr_matrix[,ncol(param_ssr_matrix)])) 
   par_str_ini = c(param_ssr_matrix[index_selection,-ncol(param_ssr_matrix)],gamma_grid[combinations[index_selection,1]],quantile(state_variable,probs = c_grid[combinations[index_selection,2]]))
   
   # Estimate by NLS LM
   nls_str_LM = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+
                             -rho*b1*x4-rho*b2*x5-rho*b3*x6
                          +(1-exp(-(gamma/sd(state_variable[2:obs]))*((s_var-c)^2)))*(b4 + b5*x1 + b6*x2 + b7*x3)
                          -rho*(1-exp(-(gamma/sd(state_variable[2:obs]))*((s_var_l-c)^2)))*(b4+ b5*x4 + b6*x5 + b7*x6)
                          , 
                          data = list(y = data_str_interest$interest_rate[2:obs], y_l = data_str_interest$interest_rate[1:obs-1], x1 = data_str_interest$interest_rate.l1[2:obs],
                                      x2 = data_str_interest$inflation_gap[2:obs], x3 = data_str_interest$output_gap[2:obs], x4 = data_str_interest$interest_rate.l1[1:obs-1], 
                                      x5 = data_str_interest$inflation_gap[1:obs-1], x6 = data_str_interest$output_gap[1:obs-1],
                                      s_var = state_variable[2:obs], s_var_l = state_variable[1:obs-1]),
                          start = list(rho = par_str_ini[1], b0 = par_str_ini[2], b1 = par_str_ini[3], b2 = par_str_ini[4], b3= par_str_ini[5], b4= par_str_ini[6], b5= par_str_ini[7], b6 = par_str_ini[8],b7 = par_str_ini[9],  
                                       gamma = par_str_ini[10], c = par_str_ini[11]),
                          
                          lower = c(0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf, -Inf,-Inf,gamma_grid[1],quantile(state_variable,probs = c_grid[1])),
                          upper = c(Inf, Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf, Inf,quantile(state_variable,probs = c_grid[length(c_grid)])),
                          control = nls.lm.control(maxiter = 200),nls.control(minFactor = 1/4096, warnOnly = TRUE), jac = NULL), silent = TRUE)
   if(try(summary(nls_str_LM),silent = TRUE)[2] %in% "try-error" || class(nls_str_LM) %in% "try-error"){next} else {
      convergence_message=if(nls_str_LM$convInfo$isConv == TRUE){convergence_message = "converged"} else {"no convergence"}
      #Thing we need to get from the model:
      # Name of the transition variable 
      name_state_variable = names_str_interest[i]
      # Parameters (linear and nonlinear part) 
      parameters_aux = nls_str_LM$m$getPars()
      #Std error of parameter estimates
      std_parameters_aux = summary(nls_str_LM)$coefficients[,2]
      #Transition Function (dated and ordered)
      f_estr1_aux = 1-exp(-(parameters_aux[length(parameters_aux)-1]/sd(state_variable[2:obs]))*((state_variable[2:obs]-parameters_aux[length(parameters_aux)])^2))
      #plot(state_variable[2:obs],f_estr1_aux, col = 'blue')
      f_estr1_aux_ordered = sort(f_estr1_aux, decreasing = FALSE)
      #plot(f_estr1_aux_ordered, col = 'red')
      #Predicted values: 
      predicted_str_aux = predict(nls_str_LM)
      #BIC (to compare against linear models) 
      BIC_str_aux = BIC(nls_str_LM)
      BIC_criteria = if(BIC_str_aux >= BIC_interest_nls){BIC_criteria = "Linear"} else {BIC_criteria = "No Linear"}
      #AIC (to compare against alternative STR especifications)
      AIC_str_aux = AIC(nls_str_LM)
      AIC_Criteria = if(AIC_str_aux >= AIC_interest_nls){AIC_Criteria = "Linear"} else {AIC_Criteria = "No Linear"}
      # Residuals of the model (with 0 at the beginning): 
      residuals_str_aux = resid(nls_str_LM)
      # SSR (necessary for evaluation stage)
      ssr_str_aux = sum((predict(nls_str_LM) -data_str_interest$interest_rate[2:obs])^2)
      SSR_criteria = if(ssr_str_aux >= ssr_1_interest_nls){SSR_criteria = "Linear"} else {SSR_criteria = "No Linear"}
      #Gradient for autocorrelation test 
      gradient_aux = nls_str_LM$m$gradient()
      #Test of restriction of nonlinear parameters = 0
      F_test_nlparam = ((ssr_1_interest_nls-ssr_str_aux)/(summary(interest_eq_nls)$df[2] -summary(nls_str_LM)$df[2]))/(ssr_str_aux/(summary(nls_str_LM)$df[2]))
      pval_f_nlparam = pf(F_test_nlparam,summary(interest_eq_nls)$df[2]-summary(nls_str_LM)$df[2],summary(nls_str_LM)$df[2],lower.tail = FALSE)
      #Test of autocorrelation of residuals 
      #1. Orthogonalize errors: 
      ortho_error_autocorr = lm(residuals_str_aux ~ 1 +gradient_aux)
      reisduals_ortho_str_aux = resid(ortho_error_autocorr)
      #2. Matrix of lagged residuales (vt):
      vt_l_residuals = cbind(reisduals_ortho_str_aux,dplyr::lag(reisduals_ortho_str_aux,1),dplyr::lag(reisduals_ortho_str_aux,2),dplyr::lag(reisduals_ortho_str_aux,3),
                             dplyr::lag(reisduals_ortho_str_aux,4),dplyr::lag(reisduals_ortho_str_aux,5),dplyr::lag(reisduals_ortho_str_aux,6),dplyr::lag(reisduals_ortho_str_aux,7),
                             dplyr::lag(reisduals_ortho_str_aux,9),dplyr::lag(reisduals_ortho_str_aux,10),dplyr::lag(reisduals_ortho_str_aux,11),dplyr::lag(reisduals_ortho_str_aux,12))
      vt_l_residuals[is.na(vt_l_residuals)] <- 0
      vt_l_residuals = vt_l_residuals[,c(-1)]
      #3. Regression result and LM tests (F and Chi) 
      auto_corr_test = lm(residuals_str_aux ~ gradient_aux + vt_l_residuals)
      ssr_alt = sum((predict(auto_corr_test)-residuals_str_aux)^2) 
      F_corr_str = ((ssr_str_aux - ssr_alt)/(summary(nls_str_LM)$df[2]-auto_corr_test$df.residual))/(ssr_alt/(auto_corr_test$df.residual))
      pval_f_corr = pf(F_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual,auto_corr_test$df.residual,lower.tail = FALSE)
      chi_corr_str = length(residuals_str_aux)*((ssr_str_aux - ssr_alt)/ssr_str_aux) 
      pval_chi_corr = pchisq(chi_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual, lower.tail = FALSE) 
      #Result of Jarque-Bera and Shapiro Test 
      jb_test_aux = jarque.bera.test(residuals_str_aux)$p.value
      shapiro_aux = shapiro.test(residuals_str_aux)$p.value
      #Ramsey Test
      ols_str_ramsey = lm(residuals_str_aux ~ gradient_aux + predicted_str_aux + predicted_str_aux^2 + predicted_str_aux^3)
      ssr_ols_str_ramsey = sum((predict(ols_str_ramsey) - residuals_str_aux)^2)
      #F test statistic:
      f_ramsey_ols = ((ssr_str_aux - ssr_ols_str_ramsey)/(summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual))/(ssr_ols_str_ramsey/(ols_str_ramsey$df.residual))
      pval_f_str_ramsey = pf(f_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,ols_str_ramsey$df.residual, lower.tail = FALSE)
      #Chi test statistic 
      chi_ramsey_ols = length(residuals_str_aux)*(ssr_str_aux - ssr_ols_str_ramsey)/ssr_str_aux
      pval_chi_str_ramsey = pchisq(chi_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,lower.tail = FALSE)
      #No additive Linearity Test - LM: 
      no_ad_nolin_results_direct = NULL
      for(k in eval(ncol(data_str_interest_model)+1):eval(ncol(data_str_interest_model)+length(selected_variables_direct))){
         state_2 =  data_str_interest[2:nrow(data_str_interest),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux + eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_direct = rbind(no_ad_nolin_results_direct,no_ad_nolin_results_aux)
      }
      
      no_ad_nolin_results_state = NULL
      for(k in eval(ncol(data_str_interest_model)+length(selected_variables_direct)+1):eval(ncol(data_str_interest_model)+ncol(data_str_interest_state))){
         state_2 =  data_str_interest[2:nrow(data_str_interest),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux +state_2+eval(state_2^2)+eval(state_2^3)+ eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_state = rbind(no_ad_nolin_results_state,no_ad_nolin_results_aux)
      }
      no_ad_nolin_results_direct = na.omit(no_ad_nolin_results_direct)
      no_ad_nolin_results_state = na.omit(no_ad_nolin_results_state)
      if(any(no_ad_nolin_results_direct[,2] <= 0.05) == TRUE){no_ad_d_nl_5_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_direct[,2] <= 0.01) == TRUE){no_ad_d_nl_1_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_1_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.05) == TRUE){no_ad_s_nl_5_pc = length(which(no_ad_nolin_results_state[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.01) == TRUE){no_ad_s_nl_1_pc = length(which(no_ad_nolin_results_state[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_1_pc="No Additive Nonlinearity"}
      #Parameter Constancy Test - LM: cannot be performed do to small sample size problem (low DF)
      #Arch test: 
      noarch_test =  ArchTest(residuals_str_aux, lags = 12)
      arch_chi_statistic = noarch_test$statistic
      arch_chi_pval = noarch_test$p.value
      
      #Final Output
      model_info = rbind(convergence_message,BIC_str_aux,BIC_criteria, AIC_str_aux, AIC_Criteria, ssr_str_aux,SSR_criteria,F_corr_str,pval_f_corr,chi_corr_str,pval_chi_corr,
                         jb_test_aux,shapiro_aux,f_ramsey_ols,pval_f_str_ramsey,chi_ramsey_ols,pval_chi_str_ramsey,no_ad_d_nl_5_pc,no_ad_d_nl_1_pc,no_ad_s_nl_5_pc,no_ad_s_nl_1_pc,arch_chi_statistic,arch_chi_pval, 
                         F_test_nlparam, pval_f_nlparam) 
      colnames(model_info) = name_state_variable
      rownames(model_info) = c("Convergence", "BIC", "BIC Criteria", "AIC", "AIC Criteria", "SSR", "SSR Criteria", "F test corr", "Pval F test", "Chi2 test corr", "Pval Chi2 test",
                               "JB test", "Shapiro test", "F test Ramsey", "Pval F Ramsey", "Chi2 test Ramsey", "Pval Chi2 Ramsey", "No Ad NL D 5%","No Ad D NL 1%","No Ad NL S 5%","No Ad S NL 1%", "ARCH test stat", "ARCH test pval", 
                               "F test nl", "Pval F nl")
      Transition_function = cbind(f_estr1_aux, f_estr1_aux_ordered, state_variable[2:obs])
      colnames(Transition_function) = c("ESTR", "ESTR sorted", "State_Variable")
      final_output = list(Model_Info = model_info,Model = nls_str_LM, Transition_function = Transition_function)
      ESTR_INTEREST[[eval(i-ncol(data_str_interest_model))]] <- final_output
   }
}
names(ESTR_INTEREST) <- names_str_interest[eval(ncol(data_str_interest_model)+1):length(names_str_interest)]
proc.time() - ptm 

#Results: of model info 
ESTR_INTEREST_model_results = NULL
for(j in 1:length(ESTR_INTEREST)){aux_info = ESTR_INTEREST[[j]]$Model_Info
ESTR_INTEREST_model_results = cbind(ESTR_INTEREST_model_results,aux_info)
}
#Selected models According to BIC/AIC criteria,  
ESTR_INTEREST_selected_models = ESTR_INTEREST_model_results[,
                                                              which(ESTR_INTEREST_model_results['Convergence',] == "converged" & 
                                                                       (ESTR_INTEREST_model_results['BIC Criteria',] == "No Linear" |  ESTR_INTEREST_model_results['AIC Criteria',] == "No Linear")
                                                                    & (ESTR_INTEREST_model_results['JB test',] >= 0.05 | ESTR_INTEREST_model_results['Shapiro test',] >= 0.05)
                                                                    &  ESTR_INTEREST_model_results['Pval F nl',] <= 0.05 &  ESTR_INTEREST_model_results['Pval F test',] >= 0.05 &  ESTR_INTEREST_model_results['Pval Chi2 test',] >= 0.05
                                                                    &  ESTR_INTEREST_model_results['Pval F Ramsey',] >= 0.05 & ESTR_INTEREST_model_results['Pval Chi2 Ramsey',] >= 0.05
                                                                    & ESTR_INTEREST_model_results['ARCH test pval',] >= 0.05 & (ESTR_INTEREST_model_results['No Ad D NL 1%',] == "No Additive Nonlinearity" | ESTR_INTEREST_model_results['No Ad D NL 1%',] < 0.25) 
                                                                    & (ESTR_INTEREST_model_results['No Ad S NL 1%',] == "No Additive Nonlinearity" | ESTR_INTEREST_model_results['No Ad S NL 1%',] < 0.25))]

write.csv(ESTR_INTEREST_selected_models, "ESTR_INTEREST_selected.csv")

   
LSTR2_INTEREST = list()
ptm <- proc.time()
#For Loop for estimates of LSTR1 function 
for(i in eval(ncol(data_str_interest_model)+1):eval(ncol(data_str_interest_model)+ncol(data_str_interest_state))){
   #Declare state variable
   state_variable = data_str_interest[,i]  
   #Grid for Gamma and C parameters
   gamma_grid = seq(0.1*sd(state_variable[2:obs]), 32*sd(state_variable[2:obs])^2, length.out = 32) #Reescaled by transition variable 
   c1_grid = seq(0.1,0.5,0.04) #c1 parameter, from 10 th percentile to 50th percentile 
   c2_grid = seq(0.55, 0.95, 0.04) #c2 parameter, from 50th percentile to 90th percentile 
   combinations = as.data.frame(crossing(gamma = 1:length(gamma_grid), c1 = 1:length(c1_grid), c2 = 1:length(c2_grid)))
   dim(combinations)
   
   #LSTR 2 estimation
   ssr_linear_str = NULL
   parameters_linear_str = NULL
   get_init_param_int <- function(x){
      #Combinations of parameters
      gamma_cond <- gamma_grid[combinations[x,1]]
      c1_cond <-quantile(state_variable[2:obs], probs = c1_grid[combinations[x,2]]) #State variable transition point
      c2_cond <-quantile(state_variable[2:obs], probs = c2_grid[combinations[x,2]]) #State variable transition point
      # Compute transition function for the state variable
      f_lstr2 = 1/(1+exp(-(gamma_cond/sd(state_variable[2:obs])^2)*(state_variable-c1_cond)*(state_variable-c2_cond))) #This becomes a variable 
      # Estimate model by OLS (without rho term) to get initial parameters for NLS estimation
      ols_aux = lm(interest_rate ~ 1 + interest_rate.l1+inflation_gap + output_gap  
                   +f_lstr2+ eval(inflation_gap*f_lstr2) + eval(output_gap*f_lstr2) + eval(interest_rate.l1*f_lstr2),
                   data = cbind(data_str_interest,f_lstr2))
      par_ols = coef(ols_aux)
      par_ols[is.na(par_ols)]<-0
      # Add variable to the estimation of linear model (by NLS):
      #From x7 to x9 are the STR interactions of model variables with state variables
      #From x11 to x12 are lagged interactions with rho term
      linear_str = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+
                                -rho*b1*x4-rho*b2*x5-rho*b3*x6
                             +b4*x7 + b5*x8 + b6*x9 + b7*x10
                             -rho*b4*x11 - rho*b5*x12 - rho*b6*x13 - rho*b7*x14
                             , 
                             start = list(rho = 0.03, b0 = 0.0, b1 = 0.9, b2 = 0.4, b3= 0.15, 
                                          b4 = par_ols[5], b5 = par_ols[6], b6 = par_ols[7], b7 = par_ols[8]),
                             data = list(y = data_str_interest$interest_rate[2:obs], y_l = data_str_interest$interest_rate[1:obs-1], x1 = data_str_interest$interest_rate.l1[2:obs],
                                         x2 = data_str_interest$inflation_gap[2:obs], x3 = data_str_interest$output_gap[2:obs], x4 = data_str_interest$interest_rate.l1[1:obs-1], 
                                         x5 = data_str_interest$inflation_gap[1:obs-1], x6 = data_str_interest$output_gap[1:obs-1],
                                         x7 = f_lstr2[2:obs], x8 = data_str_interest$interest_rate.l1[2:obs]*f_lstr2[2:obs],x9 = data_str_interest$inflation_gap[2:obs]*f_lstr2[2:obs], x10 = data_str_interest$output_gap[2:obs]*f_lstr2[2:obs],
                                         x11 = f_lstr2[1:obs-1], x12 = data_str_interest$interest_rate.l1[1:obs-1]*f_lstr2[1:obs-1],x13 = data_str_interest$inflation_gap[1:obs-1]*f_lstr2[1:obs-1], x14 = data_str_interest$output_gap[1:obs-1]*f_lstr2[1:obs-1]),
                             trace =  TRUE, jac = NULL, control = nls.lm.control(maxiter = 200), nls.control(minFactor = 1/4096, warnOnly = TRUE)), silent = TRUE)
      if(try(summary(linear_str), silent = TRUE)[2] %in% "try-error"||class(linear_str) %in% "try-error"){return(c(rep(NA, 10)))} else {
         # Get the Residual Sum of Squares of the Regression
         ssr_aux = sum((predict(linear_str) - data_str_interest$interest_rate[2:obs]) ^ 2)  ## residual sum of squares
         ssr_aux = sum((predict(linear_str) - data_str_interest$interest_rate[2:obs]) ^ 2)  ## residual sum of squares
         # Get the model parameters
         parameters_aux = linear_str$m$getPars()
         parameters_ssr_linear_str = rbind(c(), c(parameters_aux, ssr_aux))
         return(parameters_ssr_linear_str)
      }
   }
   system.time({param_ssr_matrix = foreach(i = 1:nrow(combinations), .combine = rbind, .packages = c('doParallel', 'minpack.lm')) %dopar% {get_init_param_int(i)}})
   param_ssr_matrix = na.omit(param_ssr_matrix)
   index_selection = which(param_ssr_matrix[,ncol(param_ssr_matrix)] == min(param_ssr_matrix[,ncol(param_ssr_matrix)])) 
   par_str_ini = c(param_ssr_matrix[index_selection,-ncol(param_ssr_matrix)],gamma_grid[combinations[index_selection,1]],quantile(state_variable,probs = c1_grid[combinations[index_selection,2]]), quantile(state_variable,probs = c2_grid[combinations[index_selection,3]]))
   # Estimate by NLS LM
   nls_str_LM = try(nlsLM(y ~ (1-rho)*b0+rho*y_l+b1*x1+b2*x2+b3*x3+
                             -rho*b1*x4-rho*b2*x5-rho*b3*x6
                          +(1/(1+exp(-(gamma/sd(state_variable[2:obs])^2)*(s_var-c1)*(s_var-c2))))*(b4 +b5*x1 + b6*x2 + b7*x3)
                          -rho*(1/(1+exp(-(gamma/sd(state_variable[2:obs])^2)*(s_var_l-c1)*(s_var_l-c2))))*(b4+b5*x4 + b6*x5 + b7*x6)
                          , 
                          data = list(y = data_str_interest$interest_rate[2:obs], y_l = data_str_interest$interest_rate[1:obs-1], x1 = data_str_interest$interest_rate.l1[2:obs],
                                      x2 = data_str_interest$inflation_gap[2:obs], x3 = data_str_interest$output_gap[2:obs], x4 = data_str_interest$interest_rate.l1[1:obs-1], 
                                      x5 = data_str_interest$inflation_gap[1:obs-1], x6 = data_str_interest$output_gap[1:obs-1],
                                      s_var = state_variable[2:obs], s_var_l = state_variable[1:obs-1]),
                          start = list(rho = par_str_ini[1], b0 = par_str_ini[2], b1 = par_str_ini[3], b2 = par_str_ini[4], b3= par_str_ini[5], b4= par_str_ini[6], b5= par_str_ini[7], b6 = par_str_ini[8], b7 = par_str_ini[9],
                                       gamma = par_str_ini[10], c1 = par_str_ini[11], c2 = par_str_ini[12]),
                          
                          lower = c(0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf, -Inf,-Inf,gamma_grid[1],quantile(state_variable,probs = c1_grid[1]), quantile(state_variable,probs = c2_grid[1])),
                          upper = c(Inf, Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf, Inf,quantile(state_variable,probs = c1_grid[length(c1_grid)]),quantile(state_variable,probs = c2_grid[length(c2_grid)])),
                          control = nls.lm.control(maxiter = 200),nls.control(minFactor = 1/4096, warnOnly = TRUE), jac = NULL), silent = TRUE)
   if(try(summary(nls_str_LM),silent = TRUE)[2] %in% "try-error" || class(nls_str_LM) %in% "try-error"){next} else {
      convergence_message=if(nls_str_LM$convInfo$isConv == TRUE){convergence_message = "converged"} else {"no convergence"}
      #Thing we need to get from the model:
      # Name of the transition variable 
      name_state_variable = names_str_interest[i]
      # Parameters (linear and nonlinear part) 
      parameters_aux = nls_str_LM$m$getPars()
      #Std error of parameter estimates
      std_parameters_aux = summary(nls_str_LM)$coefficients[,2]
      #Transition Function (dated and ordered)
      f_lstr2_aux = 1/(1+exp(-(parameters_aux[length(parameters_aux)-2]/sd(state_variable[2:obs])^2)*(state_variable[2:obs]-parameters_aux[length(parameters_aux)-1])*(state_variable[2:obs]-parameters_aux[length(parameters_aux)])))
      #plot(f_lstr2_aux, type = 'l', col = 'blue')
      #plot(state_variable[2:obs],f_lstr2_aux)
      f_lstr2_aux_ordered = sort(f_lstr2_aux, decreasing = FALSE)
      #plot(f_lstr2_aux_ordered, col = 'red')
      #Predicted values: 
      predicted_str_aux = predict(nls_str_LM)
      #BIC (to compare against linear models) 
      BIC_str_aux = BIC(nls_str_LM)
      BIC_criteria = if(BIC_str_aux >= BIC_interest_nls){BIC_criteria = "Linear"} else {BIC_criteria = "No Linear"}
      #AIC (to compare against alternative STR especifications)
      AIC_str_aux = AIC(nls_str_LM)
      AIC_Criteria = if(AIC_str_aux >= AIC_interest_nls){AIC_Criteria = "Linear"} else {AIC_Criteria = "No Linear"}
      # Residuals of the model (with 0 at the beginning): 
      residuals_str_aux = resid(nls_str_LM)
      # SSR (necessary for evaluation stage)
      ssr_str_aux = sum((predict(nls_str_LM) -data_str_interest$interest_rate[2:obs])^2)
      SSR_criteria = if(ssr_str_aux >= ssr_1_interest_nls){SSR_criteria = "Linear"} else {SSR_criteria = "No Linear"}
      #Gradient for autocorrelation test 
      gradient_aux = nls_str_LM$m$gradient()
      #Test of restriction of nonlinear parameters = 0
      F_test_nlparam = ((ssr_1_interest_nls-ssr_str_aux)/(summary(interest_eq_nls)$df[2] -summary(nls_str_LM)$df[2]))/(ssr_str_aux/(summary(nls_str_LM)$df[2]))
      pval_f_nlparam = pf(F_test_nlparam,summary(interest_eq_nls)$df[2]-summary(nls_str_LM)$df[2],summary(nls_str_LM)$df[2],lower.tail = FALSE)
      #Test of autocorrelation of residuals 
      #1. Orthogonalize errors: 
      ortho_error_autocorr = lm(residuals_str_aux ~ 1 +gradient_aux)
      reisduals_ortho_str_aux = resid(ortho_error_autocorr)
      #2. Matrix of lagged residuales (vt):
      vt_l_residuals = cbind(reisduals_ortho_str_aux,dplyr::lag(reisduals_ortho_str_aux,1),dplyr::lag(reisduals_ortho_str_aux,2),dplyr::lag(reisduals_ortho_str_aux,3),
                             dplyr::lag(reisduals_ortho_str_aux,4),dplyr::lag(reisduals_ortho_str_aux,5),dplyr::lag(reisduals_ortho_str_aux,6),dplyr::lag(reisduals_ortho_str_aux,7),
                             dplyr::lag(reisduals_ortho_str_aux,9),dplyr::lag(reisduals_ortho_str_aux,10),dplyr::lag(reisduals_ortho_str_aux,11),dplyr::lag(reisduals_ortho_str_aux,12))
      vt_l_residuals[is.na(vt_l_residuals)] <- 0
      vt_l_residuals = vt_l_residuals[,c(-1)]
      #3. Regression result and LM tests (F and Chi) 
      auto_corr_test = lm(residuals_str_aux ~ gradient_aux + vt_l_residuals)
      ssr_alt = sum((predict(auto_corr_test)-residuals_str_aux)^2) 
      F_corr_str = ((ssr_str_aux - ssr_alt)/(summary(nls_str_LM)$df[2]-auto_corr_test$df.residual))/(ssr_alt/(auto_corr_test$df.residual))
      pval_f_corr = pf(F_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual,auto_corr_test$df.residual,lower.tail = FALSE)
      chi_corr_str = length(residuals_str_aux)*((ssr_str_aux - ssr_alt)/ssr_str_aux) 
      pval_chi_corr = pchisq(chi_corr_str,summary(nls_str_LM)$df[2]-auto_corr_test$df.residual, lower.tail = FALSE) 
      #Result of Jarque-Bera and Shapiro Test 
      jb_test_aux = jarque.bera.test(residuals_str_aux)$p.value
      shapiro_aux = shapiro.test(residuals_str_aux)$p.value
      #Ramsey Test
      ols_str_ramsey = lm(residuals_str_aux ~ gradient_aux + predicted_str_aux + predicted_str_aux^2 + predicted_str_aux^3)
      ssr_ols_str_ramsey = sum((predict(ols_str_ramsey) - residuals_str_aux)^2)
      #F test statistic:
      f_ramsey_ols = ((ssr_str_aux - ssr_ols_str_ramsey)/(summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual))/(ssr_ols_str_ramsey/(ols_str_ramsey$df.residual))
      pval_f_str_ramsey = pf(f_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,ols_str_ramsey$df.residual, lower.tail = FALSE)
      #Chi test statistic 
      chi_ramsey_ols = length(residuals_str_aux)*(ssr_str_aux - ssr_ols_str_ramsey)/ssr_str_aux
      pval_chi_str_ramsey = pchisq(chi_ramsey_ols,summary(nls_str_LM)$df[2]-ols_str_ramsey$df.residual,lower.tail = FALSE)
      #No additive Linearity Test - LM: 
      no_ad_nolin_results_direct = NULL
      for(k in eval(ncol(data_str_interest_model)+1):eval(ncol(data_str_interest_model)+length(selected_variables_direct))){
         state_2 =  data_str_interest[2:nrow(data_str_interest),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux + eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_direct = rbind(no_ad_nolin_results_direct,no_ad_nolin_results_aux)
      }
      
      no_ad_nolin_results_state = NULL
      for(k in eval(ncol(data_str_interest_model)+length(selected_variables_direct)+1):eval(ncol(data_str_interest_model)+ncol(data_str_interest_state))){
         state_2 =  data_str_interest[2:nrow(data_str_interest),k]
         no_ad_nlin_aux = lm(residuals_str_aux ~ gradient_aux +state_2+eval(state_2^2)+eval(state_2^3)+ eval(gradient_aux*state_2) + eval(gradient_aux*state_2^2)+eval(gradient_aux*state_2^3))
         ssr_no_ad = sum((predict(no_ad_nlin_aux) - residuals_str_aux)^2) 
         f_test_no_ad = ((ssr_str_aux - ssr_no_ad)/(summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2]))/(ssr_no_ad/summary(no_ad_nlin_aux)$df[2])
         pval_f_no_ad = pf(f_test_no_ad,summary(nls_str_LM)$df[2]-summary(no_ad_nlin_aux)$df[2],summary(no_ad_nlin_aux)$df[2],lower.tail = FALSE)
         no_ad_nolin_results_aux = cbind(f_test_no_ad,pval_f_no_ad) 
         no_ad_nolin_results_state = rbind(no_ad_nolin_results_state,no_ad_nolin_results_aux)
      }
      no_ad_nolin_results_direct = na.omit(no_ad_nolin_results_direct)
      no_ad_nolin_results_state = na.omit(no_ad_nolin_results_state)
      if(any(no_ad_nolin_results_direct[,2] <= 0.05) == TRUE){no_ad_d_nl_5_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_direct[,2] <= 0.01) == TRUE){no_ad_d_nl_1_pc = length(which(no_ad_nolin_results_direct[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+1):eval(ncol(data_str_inf_model)+length(selected_variables_direct)))} else {no_ad_d_nl_1_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.05) == TRUE){no_ad_s_nl_5_pc = length(which(no_ad_nolin_results_state[,2] <= 0.05))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_5_pc="No Additive Nonlinearity"}
      if(any(no_ad_nolin_results_state[,2] <= 0.01) == TRUE){no_ad_s_nl_1_pc = length(which(no_ad_nolin_results_state[,2] <= 0.01))/length(eval(ncol(data_str_inf_model)+length(selected_variables_direct)+1):eval(ncol(data_str_inf_model)+ncol(data_str_inf_state)))} else {no_ad_s_nl_1_pc="No Additive Nonlinearity"}
      #Parameter Constancy Test - LM: cannot be performed do to small sample size problem (low DF)
      #Arch test: 
      noarch_test =  ArchTest(residuals_str_aux, lags = 12)
      arch_chi_statistic = noarch_test$statistic
      arch_chi_pval = noarch_test$p.value
      
      #Final Output
      model_info = rbind(convergence_message,BIC_str_aux,BIC_criteria, AIC_str_aux, AIC_Criteria, ssr_str_aux,SSR_criteria,F_corr_str,pval_f_corr,chi_corr_str,pval_chi_corr,
                         jb_test_aux,shapiro_aux,f_ramsey_ols,pval_f_str_ramsey,chi_ramsey_ols,pval_chi_str_ramsey,no_ad_d_nl_5_pc,no_ad_d_nl_1_pc,no_ad_s_nl_5_pc,no_ad_s_nl_1_pc,arch_chi_statistic,arch_chi_pval, 
                         F_test_nlparam, pval_f_nlparam) 
      colnames(model_info) = name_state_variable
      rownames(model_info) = c("Convergence", "BIC", "BIC Criteria", "AIC", "AIC Criteria", "SSR", "SSR Criteria", "F test corr", "Pval F test", "Chi2 test corr", "Pval Chi2 test",
                               "JB test", "Shapiro test", "F test Ramsey", "Pval F Ramsey", "Chi2 test Ramsey", "Pval Chi2 Ramsey", "No Ad NL D 5%","No Ad D NL 1%","No Ad NL S 5%","No Ad S NL 1%", "ARCH test stat", "ARCH test pval", 
                               "F test nl", "Pval F nl")
      Transition_function = cbind(f_lstr2_aux, f_lstr2_aux_ordered, state_variable[2:obs])
      colnames(Transition_function) = c("LSTR1", "LSTR1 sorted", "State_Variable")
      final_output = list(Model_Info = model_info,Model = nls_str_LM, Transition_function = Transition_function)
      LSTR2_INTEREST[[eval(i-ncol(data_str_interest_model))]] <- final_output
   }
}
names(LSTR2_INTEREST) <- names_str_interest[eval(ncol(data_str_interest_model)+1):length(names_str_interest)]
proc.time() - ptm 

#Results: of model info 
LSTR2_INTEREST_model_results = NULL
for(j in 1:length(LSTR2_INTEREST)){aux_info = LSTR2_INTEREST[[j]]$Model_Info
LSTR2_INTEREST_model_results = cbind(LSTR2_INTEREST_model_results,aux_info)
}
#Selected models According to BIC/AIC criteria,  
LSTR2_INTEREST_selected_models = LSTR2_INTEREST_model_results[,
                                                            which(LSTR2_INTEREST_model_results['Convergence',] == "converged" & 
                                                                     (LSTR2_INTEREST_model_results['BIC Criteria',] == "No Linear" |  LSTR2_INTEREST_model_results['AIC Criteria',] == "No Linear")
                                                                  & (LSTR2_INTEREST_model_results['JB test',] >= 0.05 | LSTR2_INTEREST_model_results['Shapiro test',] >= 0.05)
                                                                  &  LSTR2_INTEREST_model_results['Pval F nl',] <= 0.05 &  LSTR2_INTEREST_model_results['Pval F test',] >= 0.05 &  LSTR2_INTEREST_model_results['Pval Chi2 test',] >= 0.05
                                                                  &  LSTR2_INTEREST_model_results['Pval F Ramsey',] >= 0.05 & LSTR2_INTEREST_model_results['Pval Chi2 Ramsey',] >= 0.05
                                                                  & LSTR2_INTEREST_model_results['ARCH test pval',] >= 0.05 & (LSTR2_INTEREST_model_results['No Ad D NL 1%',] == "No Additive Nonlinearity" | LSTR2_INTEREST_model_results['No Ad D NL 1%',] < 0.25) 
                                                                  & (LSTR2_INTEREST_model_results['No Ad S NL 1%',] == "No Additive Nonlinearity" | LSTR2_INTEREST_model_results['No Ad S NL 1%',] < 0.25))]

write.csv(LSTR2_INTEREST_selected_models, "LSTR2_INTEREST_selected.csv")


# --------- EStimate of GIRF of nonlinear multivariate MMT models ------------ #
# ------------------------------------------------------------------------------ 
# Koop (1996) for multivariate models:
# 1. Get the variance covariance matrix of a MMT model residuals (K*T) ~ (3*68)
# 2. Compute Cholesky factorization of empirical resid var-cov, its inverse and orthogonalize resid using Cholesky factor
# 3. Draw randomly from the independent (K*T) resid matrix the interest rate shock (K - dimensional), and use Cholesky factor to return to dependence.   
     # Poztek uses percentile 0 to 100 in this step of the empirical distribution of interest rate equation. 
     # Here the propose approach is to use the decomposition, take the percentiles of the interes rate equation and then return to dependence. 
# 4. Select shock (percentile of interest equation distribution derived in 4) and a point in time for the model variable (from t+3.....T) (two loops)
# 5. For a given horizon N (here set at 16 periods -4 years-), randomly sample with (N+1)*R (R = repetitions/10000) from the K*T innovation matrix 
# 6. Use first N random innovations to iterate on nonlinear equations to get y(vt, wt-1), and N+1 random innovations to iterate on nonlinear model to obtain y(wt-1) 
# 7. Take the difference of y(vt, wt-1) - y(wt-1) for each N. (this is the GIRF)
# 8. Standardize GIRF by initial shock (chosen in step 4). 
# 9. Repeat step 6-8 R times to form the distribution of Monte Carlo simulation 
# 10.Compute measures of interest (i.e. median, 5th, 20th, 75th, 80th pctiles)

#Initialize and register parallel backend 
stopCluster(cl)
n_cores = detectCores()-1
registerDoParallel(cores = n_cores)
cl = makeCluster(n_cores)
registerDoSNOW(cl)

# GIRF <- function(model_inf, model_ygap, model_int, model_des = list(), N, reps_mc, reps_boot){
#    
# }

resid_inf_girf = resid(LSTR1_INF$IICV$Model)
#resid_inf_girf = resid(inflation_eq_nls)[c(-1)]
resid_outputgap_girf = resid(outputgap_eq_nls)[c(-1)]
#resid_outputgap_girf = resid(LSTR1_OUTPUTGAP$brent_gap.1$Model)
resid_interest_girf = resid(interest_eq_nls)[c(-1)]
#resid_interest_girf = resid(ESTR_INTEREST$ied_yr_primary$Model)
rmse_interest = sqrt(mean(resid_interest_girf^2))

histogram(resid_interest_girf)
ecdf(resid_interest_girf)
# 0. Definition of interest rate shock based on percentiles of empirical distribution of interest rate equation residuals: 
mp_shock_1 = quantile(resid_interest_girf, probs = seq(0, 1, 0.01)) #101 Shocks
histogram(mp_shock_1)

# 1. Var-Cov matrix of K (3) - dimensional innovation 
   #Matrix of residuals:
   residuals_matrix = cbind(resid_inf_girf, resid_outputgap_girf, resid_interest_girf)
   #Var-COV:
   emprirical_vcov = cov(residuals_matrix)
# 2. Cholesky decomposition of Var-Cov:  
   #2. 1 
   choles_decom = t(chol(emprirical_vcov))
   if(det(choles_decom) != 0){print("todo bien")} else {print("todo mal")}
   isTRUE(choles_decom%*%t(choles_decom) == emprirical_vcov)
   # 2.2 Inverse of Cholesky decomposition (factor)
   inv_chol_decom = inv(choles_decom)
   # Check inverse calculation (fails) 
   isTRUE(inv_chol_decom%*%choles_decom == diag(1, 3,3))
   #2.3 Obtain the independent matrix of innovations using the colesky factor  
   indep_residuals_matrix = t(inv_chol_decom%*%t(residuals_matrix))
# 3. Bootstrap from independent K-dimensional innovation 10000 times to have a distribution of orthogonal monetary policy shock  
sample_test_dep = matrix(NA, 10000,3)
for(i in 1:10000){
   sample_test_indep = cbind(sample(indep_residuals_matrix[,1], size = 1, replace = TRUE),sample(indep_residuals_matrix[,2], size = 1, replace = TRUE),sample(indep_residuals_matrix[,3], size = 1, replace = TRUE))
   sample_test_dep_aux = t(choles_decom%*%t(sample_test_indep))
   sample_test_dep[i,] = sample_test_dep_aux
}
histogram(sample_test_dep[,3])
mp_shock_2 = quantile(sample_test_dep[,3], probs = seq(0, 1, 0.01)) #101 policy shocks (alternative is to take all 10m reps)
mp_shock_2[11]*100
mp_shock_2[26]*100
mp_shock_2[51]*100
mp_shock_2[76]*100
mp_shock_2[91]*100
histogram(mp_shock_2)
   

# 4. Select shock and history for computation of impulse response function    
    #Base are data_str_inf_model, data_str_outputgap_model and data_str_inf_model
    data_girf_inf = data_str_inf_model[-1, -2]
    data_girf_outputgap = data_str_outputgap_model[-1,-2]
    data_girf_interest = data_str_interest_model[-1,-2]
    real_neutral_rate = data_xts$real_natural_rate_col["2003-03-01/2019-12-01"]
    
    #Parameters of each model: 
    param_inf = c(LSTR1_INF$IICV$Model$m$getPars())
    state_var_inf_girf = LSTR1_INF$IICV$Transition_function[,3]
    plot(state_var_inf_girf, type = 'l')
    #param_inf = c(inflation_eq_nls$m$getPars())
    #state_var_inf_girf = NULL
    
    param_ygap = c(outputgap_eq_nls$m$getPars())
    #param_ygap = c(LSTR1_OUTPUTGAP$brent_gap.1$Model$m$getPars())
    state_var_ygap_girf = NULL
    #state_var_ygap_girf = LSTR1_OUTPUTGAP$brent_gap.1$Transition_function[,3]
    
    param_int = c(interest_eq_nls$m$getPars())
    state_var_int_girf = NULL
    #param_int = c(ESTR_INTEREST$ied_yr_primary$Model$m$getPars())
    #state_var_int_girf = ESTR_INTEREST$ied_yr_primary$Transition_function[,3]
    
# 6.  Iterate on the nonlienar models: 
     #Shocked Series: 

   
   #GIRF_T = function(x)
   #Data for propagation:    
   #data_girf_inf = data_str_inf_model[-1, -2]
   #data_girf_outputgap = data_str_outputgap_model[-1,-2]
   #data_girf_interest = data_str_interest_model[-1,-2]
   #real_neutral_rate = data_xts$real_natural_rate_col["2003-03-01/2019-12-01"]
   
   #State variables historic data for propagation: 
   #state_var_inf_girf = LSTR1_INF$output_gap.1$Transition_function[,3]
   #state_var_int_girf = LSTR1_INTEREST$ICC$Transition_function[,3]
   
   
   #data_girf_inf = data_girf_inf[x:nrow(data_girf_inf),]
   #data_girf_outputgap = data_girf_outputgap[x:nrow(data_girf_outputgap),]
   #data_girf_interest = data_girf_interest[x:nrow(data_girf_interest),]
   #real_neutral_rate = real_neutral_rate[x:nrow(real_neutral_rate),]
   #state_var_inf_girf = state_var_inf_girf[x:length(state_var_inf_girf)]
   #state_var_int_girf = state_var_int_girf[x:length(state_var_int_girf)]
   
   
   #Given AR(1) error, shock and loops start with t=2 (1 missing first observation do to lag)
   #Code is divided in three blocks: N = 0 calculation, N = 1 to N = 2 calculation and finally N = 3 to N = 16
   
   GIRF_BOOT <- function(N, shock, model_espec = list(model_inf = c("LSTR1"), model_ygap = c("LINEAR"), model_int= c("LINEAR")), state_spec = list(type = c("Direct"), name =  c("outputgap"), lag = 0)){
      ##Parameters for testing the function
         #N = 16
         #shock = mp_shock_2[2]
         #model_espec = list(model_inf = c("LSTR1"), model_ygap = c("LINEAR"), model_int= c("ESTR"))
         #state_spec = list(type = c("Indirect"), name =  c("outputgap"), lag = 0)
      #Random sampling of monetary policy shock (orthogonalized) of size N+1 (residuals taken in triplets with no replacements) 
      sample_test_indep_N1 = cbind(sample(indep_residuals_matrix[,1], size = N, replace = TRUE),sample(indep_residuals_matrix[,2], size = N, replace = TRUE),sample(indep_residuals_matrix[,3], size = N, replace = TRUE))
      sample_test_dep_N1 = t(choles_decom%*%t(sample_test_indep_N1))
      vt_mp_shock_N = sample_test_dep_N1[2:N,3] 
      vt_mp_shock_N1 = sample_test_dep_N1[1:N,3]
      
      #Shocked Series Block
      #1. Initialize the series: 
      data_girf_inf_s = data_girf_inf
      data_girf_outputgap_s = data_girf_outputgap
      data_girf_interest_s = data_girf_interest
      real_neutral_rate_s = real_neutral_rate
      if(model_espec$model_ygap == "LINEAR"){ygap_0 = data_girf_outputgap_s[1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[2,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[1,-1]%*%param_ygap[3:7]
      } else if(model_espec$model_ygap == "LSTR1"){ygap_0 = data_girf_outputgap_s[1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[2,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[2]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_s[2,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[1]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_s[1,-1]%*%param_ygap[9:13])
      } else if(model_espec$model_ygap == "LSTR2"){ygap_0 = data_girf_outputgap_s[1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[2,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[2]-param_ygap[15])*(state_var_ygap_girf[2]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_s[2,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[1]-param_ygap[15])*(state_var_ygap_girf[1]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_s[1,-1]%*%param_ygap[9:13])
      } else if(model_espec$model_ygap == "ESTR"){ygap_0 = data_girf_outputgap_s[1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[2,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[1,-1]%*%param_ygap[3:7] + ((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[2]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_s[2,-1]%*%param_ygap[9:13]) - param_ygap[1]*((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[1]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_s[1,-1]%*%param_ygap[9:13])
      } else {print("Error espec model")}
         #replace values in data set that are affected
         data_girf_inf_s[2,3]<-ygap_0 #outputgap in inflation eq
         data_girf_inf_s[3,4]<-ygap_0 #lagged outputgap in inflation eq
         data_girf_outputgap_s[2,1]<-ygap_0 #outputgap in outputgap eq
         data_girf_outputgap_s[3,2]<-ygap_0 #lagged outputgap in outputgao eq
         data_girf_interest_s[2,3] <- ygap_0  #output gap in interest eq
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_ygap_girf[2] <- ygap_0}
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_ygap_girf[1] <- ygap_0}
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_inf_girf[2] <- ygap_0}
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_inf_girf[1] <- ygap_0}
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_int_girf[2] <- ygap_0}
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_int_girf[1] <- ygap_0}
         
      
      if(model_espec$model_inf == "LINEAR"){inf_0 = data_girf_inf_s[1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[2, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[1, -1]%*%param_inf[3:8]
      } else if(model_espec$model_inf == "LSTR1"){inf_0 = data_girf_inf_s[1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[2, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[2]-param_inf[17]))))*(param_inf[9]+data_girf_inf_s[2, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[1]-param_inf[17]))))*(param_inf[9]+data_girf_inf_s[1, -1]%*%param_inf[10:15])
      } else if(model_espec$model_inf == "LSTR2"){inf_0 = data_girf_inf_s[1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[2, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[2]-param_inf[17])*(state_var_inf_girf[2]-param_inf[18]))))*(param_inf[9]+data_girf_inf_s[2, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[1]-param_inf[17])*(state_var_inf_girf[1]-param_inf[18]))))*(param_inf[9]+data_girf_inf_s[1, -1]%*%param_inf[10:15])
      } else if(model_espec$model_inf == "ESTR"){inf_0 = data_girf_inf_s[1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[2, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[1, -1]%*%param_inf[3:8] +((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[2]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_s[2, -1]%*%param_inf[10:15])-param_inf[1]*((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[1]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_s[1, -1]%*%param_inf[10:15])
      } else {print("Error espec model")}
         #Replace values in data set that are affected
         data_girf_inf_s[2,1]<-inf_0 #Inflation in inflation eq
         data_girf_inf_s[3,1]<-inf_0 #Lagged inflation in inflation eq
         data_girf_interest_s[2,2] <- inf_0 #Inflation in interest rate
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_ygap_girf[2] <- inf_0}
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_ygap_girf[1] <- inf_0}
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_inf_girf[2] <- inf_0}
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_inf_girf[1] <- inf_0}
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_int_girf[2] <- inf_0}
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_int_girf[1] <- inf_0}
         
      
      if(model_espec$model_int == "LINEAR"){int_0 = data_girf_interest_s[1,1]*param_int[1] + (1-param_int[1])*param_int[2] + data_girf_interest_s[2,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[1,c(4,2,3)]%*%param_int[3:5]+shock
      } else if(model_espec$model_int == "LSTR1"){int_0 = data_girf_interest_s[1,1]*param_int[1] + (1-param_int[1])*param_int[2] + data_girf_interest_s[2,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[2]-param_int[11]))))*(param_int[6]+data_girf_interest_s[2,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[1]-param_int[11]))))*(param_int[6]+data_girf_interest_s[1,c(4,2,3)]%*%param_int[7:9])+shock
      } else if(model_espec$model_int == "LSTR2"){int_0 = data_girf_interest_s[1,1]*param_int[1] + (1-param_int[1])*param_int[2] + data_girf_interest_s[2,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[2]-param_int[11])*(state_var_int_girf[2]-param_int[12]))))*(param_int[6]+data_girf_interest_s[2,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[1]-param_int[11])*(state_var_int_girf[1]-param_int[12]))))*(param_int[6]+data_girf_interest_s[1,c(4,2,3)]%*%param_int[7:9])+shock
      } else if(model_espec$model_int == "ESTR"){int_0 = data_girf_interest_s[1,1]*param_int[1] + (1-param_int[1])*param_int[2] + data_girf_interest_s[2,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[1,c(4,2,3)]%*%param_int[3:5] +((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[2]-param_int[11])^2)))*(param_int[6]+data_girf_interest_s[2,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[1]-param_int[11])^2)))*(param_int[6]+data_girf_interest_s[1,c(4,2,3)]%*%param_int[7:9])+shock
      } else {print("Error espec model")}
         #Replace variables that are affected (as the real gap is llaged in outputgap and inflation eq especs, they are not included here)
         r_gap_0 = ((1+int_0)/(1+as.numeric(data_girf_inf_s[1,1]))-1)-as.numeric(real_neutral_rate_s[2]) #Real interest rate
         data_girf_interest_s[2,1]<- int_0 #Interest rate in interest rate eq.
         data_girf_interest_s[3,4]<- int_0 #Lagged interest rate in interest rate eq.
      #Shocked Series
      inf_girf_s = inf_0
      ygap_girf_s = ygap_0
      int_girf_s = int_0
      r_gap_s = r_gap_0
      
      #Loop up to N = 2 (interest rate exerts no effect on output and and inflation)
      for(i in 3:3){
         if(model_espec$model_ygap == "LINEAR"){ygap_aux = data_girf_outputgap_s[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[i-1,-1]%*%param_ygap[3:7]
         } else if(model_espec$model_ygap == "LSTR1"){ygap_aux = data_girf_outputgap_s[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[i-1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_s[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i-1]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_s[i-1,-1]%*%param_ygap[9:13])
         } else if(model_espec$model_ygap == "LSTR2"){ygap_aux = data_girf_outputgap_s[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[i-1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[i]-param_ygap[15])*(state_var_ygap_girf[i]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_s[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[i-1]-param_ygap[15])*(state_var_ygap_girf[i-1]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_s[i-1,-1]%*%param_ygap[9:13])
         } else if(model_espec$model_ygap == "ESTR"){ygap_aux = data_girf_outputgap_s[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[i-1,-1]%*%param_ygap[3:7] + ((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_s[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i-1]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_s[i-1,-1]%*%param_ygap[9:13])
         } else {print("Error espec model")}   
            #Replace variables that are affected 
            data_girf_inf_s[i,3]<-ygap_aux #output gap in inflation eq. 
            data_girf_inf_s[i+1,4]<-ygap_aux #Lagged outputgap in inflation eq. 
            data_girf_outputgap_s[i,1]<-ygap_aux #Output gap in outputgap eq.  
            data_girf_outputgap_s[i+1,2]<-ygap_aux  #Lagged outputgap in outputgap eq. 
            data_girf_interest_s[i,3] <- ygap_aux #Outputgap in interest rate eq. 
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_ygap_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_ygap_girf[i-1] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_inf_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_inf_girf[i-1] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_int_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_int_girf[i-1] <- ygap_aux}
         
         if(model_espec$model_inf == "LINEAR"){inf_aux = data_girf_inf_s[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[i-1, -1]%*%param_inf[3:8]
         } else if(model_espec$model_inf == "LSTR1"){inf_aux = data_girf_inf_s[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[i-1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i]-param_inf[17]))))*(param_inf[9]+data_girf_inf_s[i, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i-1]-param_inf[17]))))*(param_inf[9]+data_girf_inf_s[i-1, -1]%*%param_inf[10:15])
         } else if(model_espec$model_inf == "LSTR2"){inf_aux = data_girf_inf_s[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[i-1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[i]-param_inf[17])*(state_var_inf_girf[i]-param_inf[18]))))*(param_inf[9]+data_girf_inf_s[i, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[i-1]-param_inf[17])*(state_var_inf_girf[i-1]-param_inf[18]))))*(param_inf[9]+data_girf_inf_s[i-1, -1]%*%param_inf[10:15])
         } else if(model_espec$model_inf == "ESTR"){inf_aux = data_girf_inf_s[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[i-1, -1]%*%param_inf[3:8] +((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_s[i, -1]%*%param_inf[10:15])-param_inf[1]*((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i-1]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_s[i-1, -1]%*%param_inf[10:15])
         } else {print("Error espec model")}
            #Replace variables that are affected
            data_girf_inf_s[i,1]<-inf_aux #Inflation in inflation eq.
            data_girf_inf_s[i+1,1]<-inf_aux #Lagged inflation in inflation eq. 
            data_girf_interest_s[i,2] <- inf_aux #Inflation in interest rate eq. 
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_ygap_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_ygap_girf[i-1] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_inf_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_inf_girf[i-1] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_int_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_int_girf[i-1] <- inf_aux}
      
         if(model_espec$model_int == "LINEAR"){int_aux = data_girf_interest_s[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_s[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[3:5]+vt_mp_shock_N[i-2] # vt_mp_shock_N has to be equal to 1 for starting i value
         } else if(model_espec$model_int == "LSTR1"){int_aux = data_girf_interest_s[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_s[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i]-param_int[11]))))*(param_int[6]+data_girf_interest_s[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i-1]-param_int[11]))))*(param_int[6]+data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N[i-2] # vt_mp_shock_N has to be equal to 1 for starting i value 
         } else if(model_espec$model_int == "LSTR2"){int_aux = data_girf_interest_s[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_s[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[i]-param_int[11])*(state_var_int_girf[i]-param_int[12]))))*(param_int[6]+data_girf_interest_s[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[i-1]-param_int[11])*(state_var_int_girf[i-1]-param_int[12]))))*(param_int[6]+data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N[i-2] # vt_mp_shock_N has to be equal to 1 for starting i value 
         } else if(model_espec$model_int == "ESTR"){int_aux = data_girf_interest_s[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_s[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[3:5] +((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i]-param_int[11])^2)))*(param_int[6]+data_girf_interest_s[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i-1]-param_int[11])^2)))*(param_int[6]+data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N[i-2] # vt_mp_shock_N has to be equal to 1 for starting i value 
         } else {print("Error espec model")} 
            #Replace variables that are affected 
            r_gap_aux = ((1+int_aux)/(1+as.numeric(data_girf_inf_s[i-1,1]))-1)-as.numeric(real_neutral_rate_s[i]) #Real interest rate
            data_girf_interest_s[i,1]<- int_aux #Interest rate in interest rate eq. 
            data_girf_interest_s[i+1,4]<- int_aux #Lagged interest rate in interest rate eq. 
         inf_girf_s = rbind(inf_girf_s,inf_aux)
         ygap_girf_s = rbind(ygap_girf_s, ygap_aux)
         int_girf_s = rbind(int_girf_s, int_aux)
         r_gap_s = rbind(r_gap_s,r_gap_aux)
      }
      data_girf_inf_s[4,5] <- r_gap_s[1]
      data_girf_outputgap_s[4,3] <-r_gap_s[1]
      
      #Loop for N = 3 to N = 16
      for(i in 4:eval(N+1)){
         if(model_espec$model_ygap == "LINEAR"){ygap_aux = data_girf_outputgap_s[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[i-1,-1]%*%param_ygap[3:7]
         } else if(model_espec$model_ygap == "LSTR1"){ygap_aux = data_girf_outputgap_s[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[i-1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_s[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i-1]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_s[i-1,-1]%*%param_ygap[9:13])   
         } else if(model_espec$model_ygap == "LSTR2"){ygap_aux = data_girf_outputgap_s[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[i-1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[i]-param_ygap[15])*(state_var_ygap_girf[i]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_s[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[i-1]-param_ygap[15])*(state_var_ygap_girf[i-1]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_s[i-1,-1]%*%param_ygap[9:13]) 
         } else if(model_espec$model_ygap == "ESTR"){ygap_aux = data_girf_outputgap_s[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_s[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_s[i-1,-1]%*%param_ygap[3:7] + ((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_s[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i-1]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_s[i-1,-1]%*%param_ygap[9:13])
         } else {print("Error espec model")}
            #Replace variables that are affected 
            data_girf_inf_s[i,3]<-ygap_aux #output gap in inflation eq. 
            data_girf_inf_s[i+1,4]<-ygap_aux #Lagged outputgap in inflation eq. 
            data_girf_outputgap_s[i,1]<-ygap_aux #Output gap in outputgap eq.  
            data_girf_outputgap_s[i+1,2]<-ygap_aux  #Lagged outputgap in outputgap eq. 
            data_girf_interest_s[i,3] <- ygap_aux #Outputgap in interest rate eq. 
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_ygap_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_ygap_girf[i-1] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_inf_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_inf_girf[i-1] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_int_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_int_girf[i-1] <- ygap_aux}
         
         if(model_espec$model_inf == "LINEAR"){inf_aux = data_girf_inf_s[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[i-1, -1]%*%param_inf[3:8]
         } else if(model_espec$model_inf == "LSTR1"){inf_aux = data_girf_inf_s[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[i-1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i]-param_inf[17]))))*(param_inf[9]+data_girf_inf_s[i, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i-1]-param_inf[17]))))*(param_inf[9]+data_girf_inf_s[i-1, -1]%*%param_inf[10:15])
         } else if(model_espec$model_inf == "LSTR2"){inf_aux = data_girf_inf_s[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[i-1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[i]-param_inf[17])*(state_var_inf_girf[i]-param_inf[18]))))*(param_inf[9]+data_girf_inf_s[i, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[i-1]-param_inf[17])*(state_var_inf_girf[i-1]-param_inf[18]))))*(param_inf[9]+data_girf_inf_s[i-1, -1]%*%param_inf[10:15])
         } else if(model_espec$model_inf == "ESTR"){inf_aux = data_girf_inf_s[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_s[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_s[i-1, -1]%*%param_inf[3:8] +((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_s[i, -1]%*%param_inf[10:15])-param_inf[1]*((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i-1]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_s[i-1, -1]%*%param_inf[10:15])
         } else {print("Error espec model")}   
            #Replace variables that are affected
            data_girf_inf_s[i,1]<-inf_aux #Inflation in inflation eq.
            data_girf_inf_s[i+1,1]<-inf_aux #Lagged inflation in inflation eq. 
            data_girf_interest_s[i,2] <- inf_aux #Inflation in interest rate eq. 
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_ygap_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_ygap_girf[i-1] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_inf_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_inf_girf[i-1] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_int_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_int_girf[i-1] <- inf_aux}
         
         
         if(model_espec$model_int == "LINEAR"){int_aux = data_girf_interest_s[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_s[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[3:5]+vt_mp_shock_N[i-2] # vt_mp_shock_N has to be equal to 1 for starting i value
         } else if(model_espec$model_int == "LSTR1"){int_aux = data_girf_interest_s[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_s[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i]-param_int[11]))))*(param_int[6]+data_girf_interest_s[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i-1]-param_int[11]))))*(param_int[6]+data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N[i-2] # vt_mp_shock_N has to be equal to 1 for starting i value 
         } else if(model_espec$model_int == "LSTR2"){int_aux = data_girf_interest_s[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_s[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[i]-param_int[11])*(state_var_int_girf[i]-param_int[12]))))*(param_int[6]+data_girf_interest_s[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[i-1]-param_int[11])*(state_var_int_girf[i-1]-param_int[12]))))*(param_int[6]+data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N[i-2] # vt_mp_shock_N has to be equal to 1 for starting i value 
         } else if(model_espec$model_int == "ESTR"){int_aux = data_girf_interest_s[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_s[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[3:5] +((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i]-param_int[11])^2)))*(param_int[6]+data_girf_interest_s[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i-1]-param_int[11])^2)))*(param_int[6]+data_girf_interest_s[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N[i-2] # vt_mp_shock_N has to be equal to 1 for starting i value 
         } else {print("Error espec model")}  
            #Replace variables that are affected 
            r_gap_aux = ((1+int_aux)/(1+as.numeric(data_girf_inf_s[i-1,1]))-1)-as.numeric(real_neutral_rate_s[i]) #Real interest rate. 
            data_girf_interest_s[i,1]<- int_aux #Interest rate in interest rate eq. 
            data_girf_interest_s[i+1,4]<- int_aux #Lagged interest rate in in interest rate eq. 
            data_girf_inf_s[i+1,5] <- r_gap_s[i-2] #Real interest rate in inflation equation. 
            data_girf_outputgap_s[i+1,3] <-r_gap_s[i-2] #Real interest rate in output gap equation. 
         
         inf_girf_s = rbind(inf_girf_s,inf_aux)
         ygap_girf_s = rbind(ygap_girf_s, ygap_aux)
         int_girf_s = rbind(int_girf_s, int_aux)
         r_gap_s = rbind(r_gap_s,r_gap_aux)
      }
      
      #Non Shocked Series
      #Reset propagation data of models
      data_girf_inf_ns = data_girf_inf
      data_girf_outputgap_ns = data_girf_outputgap
      data_girf_interest_ns = data_girf_interest
      real_neutral_rate_ns = real_neutral_rate
      
      if(model_espec$model_ygap == "LINEAR"){ygap_0 = data_girf_outputgap_ns[1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[2,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[1,-1]%*%param_ygap[3:7]
      } else if(model_espec$model_ygap == "LSTR1"){ygap_0 = data_girf_outputgap_ns[1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[2,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[2]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_ns[2,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[1]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_ns[1,-1]%*%param_ygap[9:13])
      } else if(model_espec$model_ygap == "LSTR2"){ygap_0 = data_girf_outputgap_ns[1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[2,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[2]-param_ygap[15])*(state_var_ygap_girf[2]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_ns[2,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[1]-param_ygap[15])*(state_var_ygap_girf[1]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_ns[1,-1]%*%param_ygap[9:13])
      } else if(model_espec$model_ygap == "ESTR"){ygap_0 = data_girf_outputgap_ns[1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[2,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[1,-1]%*%param_ygap[3:7] + ((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[2]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_ns[2,-1]%*%param_ygap[9:13]) - param_ygap[1]*((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[1]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_ns[1,-1]%*%param_ygap[9:13])
      } else {print("Error espec model")}
         #replace values in data set that are affected
         data_girf_inf_ns[2,3]<-ygap_0 #outputgap in inflation eq
         data_girf_inf_ns[3,4]<-ygap_0 #lagged outputgap in inflation eq
         data_girf_outputgap_ns[2,1]<-ygap_0 #outputgap in outputgap eq
         data_girf_outputgap_ns[3,2]<-ygap_0 #lagged outputgap in outputgao eq
         data_girf_interest_ns[2,3] <- ygap_0  #output gap in interest eq
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_ygap_girf[2] <- ygap_0}
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_ygap_girf[1] <- ygap_0}
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_inf_girf[2] <- ygap_0}
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_inf_girf[1] <- ygap_0}
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_int_girf[2] <- ygap_0}
         if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_int_girf[1] <- ygap_0}
      
      if(model_espec$model_inf == "LINEAR"){inf_0 = data_girf_inf_ns[1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[2, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[1, -1]%*%param_inf[3:8]
      } else if(model_espec$model_inf == "LSTR1"){inf_0 = data_girf_inf_ns[1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[2, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[2]-param_inf[17]))))*(param_inf[9]+data_girf_inf_ns[2, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[1]-param_inf[17]))))*(param_inf[9]+data_girf_inf_ns[1, -1]%*%param_inf[10:15])
      } else if(model_espec$model_inf == "LSTR2"){inf_0 = data_girf_inf_ns[1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[2, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[2]-param_inf[17])*(state_var_inf_girf[2]-param_inf[18]))))*(param_inf[9]+data_girf_inf_ns[2, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[1]-param_inf[17])*(state_var_inf_girf[1]-param_inf[18]))))*(param_inf[9]+data_girf_inf_ns[1, -1]%*%param_inf[10:15])
      } else if(model_espec$model_inf == "ESTR"){inf_0 = data_girf_inf_ns[1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[2, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[1, -1]%*%param_inf[3:8] +((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[2]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_ns[2, -1]%*%param_inf[10:15])-param_inf[1]*((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[1]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_ns[1, -1]%*%param_inf[10:15])
      } else {print("Error espec model")}
         #Replace values in data set that are affected
         data_girf_inf_ns[2,1]<-inf_0 #Inflation in inflation eq
         data_girf_inf_ns[3,1]<-inf_0 #Lagged inflation in inflation eq
         data_girf_interest_ns[2,2] <- inf_0 #Inflation in interest rate
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_ygap_girf[2] <- inf_0}
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_ygap_girf[1] <- inf_0}
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_inf_girf[2] <- inf_0}
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_inf_girf[1] <- inf_0}
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_int_girf[2] <- inf_0}
         if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_int_girf[1] <- inf_0}
      
      
      if(model_espec$model_int == "LINEAR"){int_0 = data_girf_interest_ns[1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[2,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[1,c(4,2,3)]%*%param_int[3:5]+0
      } else if(model_espec$model_int == "LSTR1"){int_0 = data_girf_interest_ns[1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[2,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[2]-param_int[11]))))*(param_int[6]+data_girf_interest_ns[2,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[1]-param_int[11]))))*(param_int[6]+data_girf_interest_ns[1,c(4,2,3)]%*%param_int[7:9])+0 #vt_mp_shock_N1[1]
      } else if(model_espec$model_int == "LSTR2"){int_0 = data_girf_interest_ns[1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[2,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[2]-param_int[11])*(state_var_int_girf[2]-param_int[12]))))*(param_int[6]+data_girf_interest_ns[2,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[1]-param_int[11])*(state_var_int_girf[1]-param_int[12]))))*(param_int[6]+data_girf_interest_ns[1,c(4,2,3)]%*%param_int[7:9])+0#vt_mp_shock_N1[1]
      } else if(model_espec$model_int == "ESTR"){int_0 = data_girf_interest_ns[1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[2,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[1,c(4,2,3)]%*%param_int[3:5] +((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[2]-param_int[11])^2)))*(param_int[6]+data_girf_interest_ns[2,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[1]-param_int[11])^2)))*(param_int[6]+data_girf_interest_ns[1,c(4,2,3)]%*%param_int[7:9])+0 #vt_mp_shock_N1[1]
      } else {print("Error espec model")}
         #Replace variables that are affected (as the real gap is llaged in outputgap and inflation eq especs, they are not included here)
         r_gap_0 = ((1+int_0)/(1+as.numeric(data_girf_inf_ns[1,1]))-1)-as.numeric(real_neutral_rate_ns[2]) #Real interest rate
         data_girf_interest_ns[2,1]<- int_0 #Interest rate in interest rate eq.
         data_girf_interest_ns[3,4]<- int_0 #Lagged interest rate in interest rate eq.
      
      inf_girf_ns = inf_0
      ygap_girf_ns = ygap_0
      int_girf_ns = int_0
      r_gap_ns = r_gap_0
      
      #Loop up to N = 2 (interest rate exerts no effect on output and and inflation)
      for(i in 3:3){
         if(model_espec$model_ygap == "LINEAR"){ygap_aux = data_girf_outputgap_ns[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[i-1,-1]%*%param_ygap[3:7]
         } else if(model_espec$model_ygap == "LSTR1"){ygap_aux = data_girf_outputgap_ns[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[i-1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_ns[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i-1]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_ns[i-1,-1]%*%param_ygap[9:13])
         } else if(model_espec$model_ygap == "LSTR2"){ygap_aux = data_girf_outputgap_ns[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[i-1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[i]-param_ygap[15])*(state_var_ygap_girf[i]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_ns[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[i-1]-param_ygap[15])*(state_var_ygap_girf[i-1]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_ns[i-1,-1]%*%param_ygap[9:13])
         } else if(model_espec$model_ygap == "ESTR"){ygap_aux = data_girf_outputgap_ns[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[i-1,-1]%*%param_ygap[3:7] + ((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_ns[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i-1]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_ns[i-1,-1]%*%param_ygap[9:13])
         } else {print("Error espec model")}
            #Replace variables that are affected 
            data_girf_inf_ns[i,3]<-ygap_aux #output gap in inflation eq. 
            data_girf_inf_ns[i+1,4]<-ygap_aux #Lagged outputgap in inflation eq. 
            data_girf_outputgap_ns[i,1]<-ygap_aux #Output gap in outputgap eq.  
            data_girf_outputgap_ns[i+1,2]<-ygap_aux  #Lagged outputgap in outputgap eq. 
            data_girf_interest_ns[i,3] <- ygap_aux #Outputgap in interest rate eq. 
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_ygap_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_ygap_girf[i-1] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_inf_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_inf_girf[i-1] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_int_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_int_girf[i-1] <- ygap_aux}
         
         
         if(model_espec$model_inf == "LINEAR"){inf_aux = data_girf_inf_ns[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[i-1, -1]%*%param_inf[3:8]
         } else if(model_espec$model_inf == "LSTR1"){inf_aux = data_girf_inf_ns[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[i-1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i]-param_inf[17]))))*(param_inf[9]+data_girf_inf_ns[i, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i-1]-param_inf[17]))))*(param_inf[9]+data_girf_inf_ns[i-1, -1]%*%param_inf[10:15])
         } else if(model_espec$model_inf == "LSTR2"){inf_aux = data_girf_inf_ns[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[i-1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[i]-param_inf[17])*(state_var_inf_girf[i]-param_inf[18]))))*(param_inf[9]+data_girf_inf_ns[i, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[i-1]-param_inf[17])*(state_var_inf_girf[i-1]-param_inf[18]))))*(param_inf[9]+data_girf_inf_ns[i-1, -1]%*%param_inf[10:15])   
         } else if(model_espec$model_inf == "ESTR"){inf_aux = data_girf_inf_ns[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[i-1, -1]%*%param_inf[3:8] +((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_ns[i, -1]%*%param_inf[10:15])-param_inf[1]*((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i-1]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_ns[i-1, -1]%*%param_inf[10:15])
         } else {print("Error espec model")}
            #Replace variables that are affected
            data_girf_inf_ns[i,1]<-inf_aux #Inflation in inflation eq.
            data_girf_inf_ns[i+1,1]<-inf_aux #Lagged inflation in inflation eq. 
            data_girf_interest_ns[i,2] <- inf_aux #Inflation in interest rate eq. 
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_ygap_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_ygap_girf[i-1] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_inf_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_inf_girf[i-1] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_int_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_int_girf[i-1] <- inf_aux}
         
         if(model_espec$model_int == "LINEAR"){int_aux = data_girf_interest_ns[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[3:5]+vt_mp_shock_N1[i-1] # vt_mp_shock_N has to be equal to 1 for starting i value 
         } else if(model_espec$model_int == "LSTR1"){int_aux = data_girf_interest_ns[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i]-param_int[11]))))*(param_int[6]+data_girf_interest_ns[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i-1]-param_int[11]))))*(param_int[6]+data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N1[i-1]
         } else if(model_espec$model_int == "LSTR2"){int_aux = data_girf_interest_ns[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[i]-param_int[11])*(state_var_int_girf[i]-param_int[12]))))*(param_int[6]+data_girf_interest_ns[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[i-1]-param_int[11])*(state_var_int_girf[i-1]-param_int[12]))))*(param_int[6]+data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N1[i-1]
         } else if(model_espec$model_int == "ESTR"){int_aux = data_girf_interest_ns[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[3:5] +((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i]-param_int[11])^2)))*(param_int[6]+data_girf_interest_ns[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i-1]-param_int[11])^2)))*(param_int[6]+data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N1[i-1]
         } else {print("Error espec model")}
            #Replace variables that are affected 
            r_gap_aux = ((1+int_aux)/(1+as.numeric(data_girf_inf_ns[i-1,1]))-1)-as.numeric(real_neutral_rate_ns[i]) #Real interest rate
            data_girf_interest_ns[i,1]<- int_aux #Interest rate in interest rate eq. 
            data_girf_interest_ns[i+1,4]<- int_aux #Lagged interest rate in interest rate eq. 
         
         inf_girf_ns = rbind(inf_girf_ns,inf_aux)
         ygap_girf_ns = rbind(ygap_girf_ns, ygap_aux)
         int_girf_ns = rbind(int_girf_ns, int_aux)
         r_gap_ns = rbind(r_gap_ns,r_gap_aux)
      }
      
      data_girf_inf_ns[4,5] <- r_gap_ns[1]
      data_girf_outputgap_ns[4,3] <-r_gap_ns[1]
      
      #Loop for N = 3 to N = 16
      for(i in 4:(N+1)){
         if(model_espec$model_ygap == "LINEAR"){ygap_aux = data_girf_outputgap_ns[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[i-1,-1]%*%param_ygap[3:7]
         } else if(model_espec$model_ygap == "LSTR1"){ygap_aux = data_girf_outputgap_ns[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[i-1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_ns[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i-1]-param_ygap[15]))))*(param_ygap[8]+data_girf_outputgap_ns[i-1,-1]%*%param_ygap[9:13])
         } else if(model_espec$model_ygap == "LSTR2"){ygap_aux = data_girf_outputgap_ns[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[i-1,-1]%*%param_ygap[3:7] + (1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[i]-param_ygap[15])*(state_var_ygap_girf[i]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_ns[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*(1/(1+exp(-(param_ygap[14]/sd(state_var_ygap_girf)^2)*(state_var_ygap_girf[i-1]-param_ygap[15])*(state_var_ygap_girf[i-1]-param_ygap[16]))))*(param_ygap[8]+data_girf_outputgap_ns[i-1,-1]%*%param_ygap[9:13])
         } else if(model_espec$model_ygap == "ESTR"){ygap_aux = data_girf_outputgap_ns[i-1,1]%*%param_ygap[1] + (1-param_ygap[1])*param_ygap[2]+data_girf_outputgap_ns[i,-1]%*%param_ygap[3:7] -param_ygap[1]*data_girf_outputgap_ns[i-1,-1]%*%param_ygap[3:7] + ((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_ns[i,-1]%*%param_ygap[9:13]) - param_ygap[1]*((1-exp(-(param_ygap[14]/sd(state_var_ygap_girf))*(state_var_ygap_girf[i-1]-param_ygap[15])^2)))*(param_ygap[8]+data_girf_outputgap_ns[i-1,-1]%*%param_ygap[9:13])
         } else {print("Error espec model")} 
            #Replace variables that are affected 
            data_girf_inf_ns[i,3]<-ygap_aux #output gap in inflation eq. 
            data_girf_inf_ns[i+1,4]<-ygap_aux #Lagged outputgap in inflation eq. 
            data_girf_outputgap_ns[i,1]<-ygap_aux #Output gap in outputgap eq.  
            data_girf_outputgap_ns[i+1,2]<-ygap_aux  #Lagged outputgap in outputgap eq. 
            data_girf_interest_ns[i,3] <- ygap_aux #Outputgap in interest rate eq. 
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_ygap_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_ygap_girf[i-1] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_inf_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_inf_girf[i-1] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 0){state_var_int_girf[i] <- ygap_aux}
            if(state_spec$type == "Direct" & state_spec$name == "outputgap" & state_spec$lag == 1){state_var_int_girf[i-1] <- ygap_aux}
         
         
         if(model_espec$model_inf == "LINEAR"){inf_aux = data_girf_inf_ns[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[i-1, -1]%*%param_inf[3:8]
         } else if(model_espec$model_inf == "LSTR1"){inf_aux = data_girf_inf_ns[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[i-1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i]-param_inf[17]))))*(param_inf[9]+data_girf_inf_ns[i, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i-1]-param_inf[17]))))*(param_inf[9]+data_girf_inf_ns[i-1, -1]%*%param_inf[10:15])
         } else if(model_espec$model_inf == "LSTR2"){inf_aux = data_girf_inf_ns[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[i-1, -1]%*%param_inf[3:8] +(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[i]-param_inf[17])*(state_var_inf_girf[i]-param_inf[18]))))*(param_inf[9]+data_girf_inf_ns[i, -1]%*%param_inf[10:15])-param_inf[1]*(1/(1+exp(-(param_inf[16]/sd(state_var_inf_girf)^2)*(state_var_inf_girf[i-1]-param_inf[17])*(state_var_inf_girf[i-1]-param_inf[18]))))*(param_inf[9]+data_girf_inf_ns[i-1, -1]%*%param_inf[10:15])
         } else if(model_espec$model_inf == "ESTR"){inf_aux = data_girf_inf_ns[i-1,1]%*%param_inf[1]+(1-param_inf[1])*param_inf[2]+data_girf_inf_ns[i, -1]%*%param_inf[3:8]-param_inf[1]*data_girf_inf_ns[i-1, -1]%*%param_inf[3:8] +((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_ns[i, -1]%*%param_inf[10:15])-param_inf[1]*((1-exp(-(param_inf[16]/sd(state_var_inf_girf))*(state_var_inf_girf[i-1]-param_inf[17])^2)))*(param_inf[9]+data_girf_inf_ns[i-1, -1]%*%param_inf[10:15])
         } else {print("Error espec model")}
            #Replace variables that are affected
            data_girf_inf_ns[i,1]<-inf_aux #Inflation in inflation eq.
            data_girf_inf_ns[i+1,1]<-inf_aux #Lagged inflation in inflation eq. 
            data_girf_interest_ns[i,2] <- inf_aux #Inflation in interest rate eq. 
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_ygap_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_ygap_girf[i-1] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_inf_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_inf_girf[i-1] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 0){state_var_int_girf[i] <- inf_aux}
            if(state_spec$type == "Direct" & state_spec$name == "inflation" & state_spec$lag == 1){state_var_int_girf[i-1] <- inf_aux}
            
         
         if(model_espec$model_int == "LINEAR"){int_aux = data_girf_interest_ns[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[3:5]+vt_mp_shock_N1[i-1] # vt_mp_shock_N has to be equal to 1 for starting i value 
         } else if(model_espec$model_int == "LSTR1"){int_aux = data_girf_interest_ns[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i]-param_int[11]))))*(param_int[6]+data_girf_interest_ns[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i-1]-param_int[11]))))*(param_int[6]+data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N1[i-1]
         } else if(model_espec$model_int == "LSTR2"){int_aux = data_girf_interest_ns[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[3:5] +(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[i]-param_int[11])*(state_var_int_girf[i]-param_int[12]))))*(param_int[6]+data_girf_interest_ns[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*(1/(1+exp(-(param_int[10]/sd(state_var_int_girf)^2)*(state_var_int_girf[i-1]-param_int[11])*(state_var_int_girf[i-1]-param_int[12]))))*(param_int[6]+data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N1[i-1]
         } else if(model_espec$model_int == "ESTR"){int_aux = data_girf_interest_ns[i-1,1]*param_int[1] + (1-param_int[1])*param_int[2] +data_girf_interest_ns[i,c(4,2,3)]%*%param_int[3:5]-param_int[1]*data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[3:5] +((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i]-param_int[11])^2)))*(param_int[6]+data_girf_interest_ns[i,c(4,2,3)]%*%param_int[7:9]) -param_int[1]*((1-exp(-(param_int[10]/sd(state_var_int_girf))*(state_var_int_girf[i-1]-param_int[11])^2)))*(param_int[6]+data_girf_interest_ns[i-1,c(4,2,3)]%*%param_int[7:9])+vt_mp_shock_N1[i-1]
         } else {print("Error espec model")}
            #Replace variables that are affected 
            r_gap_aux = ((1+int_aux)/(1+as.numeric(data_girf_inf_ns[i-1,1]))-1)-as.numeric(real_neutral_rate_ns[i]) #Real interest rate. 
            data_girf_interest_ns[i,1]<- int_aux #Interest rate in interest rate eq. 
            data_girf_interest_ns[i+1,4]<- int_aux #Lagged interest rate in in interest rate eq. 
            data_girf_inf_ns[i+1,5] <- r_gap_ns[i-2] #Real interest rate in inflation equation. 
            data_girf_outputgap_ns[i+1,3] <-r_gap_ns[i-2] #Real interest rate in output gap equation.
         
         inf_girf_ns = rbind(inf_girf_ns,inf_aux)
         ygap_girf_ns = rbind(ygap_girf_ns, ygap_aux)
         int_girf_ns = rbind(int_girf_ns, int_aux)
         r_gap_ns = rbind(r_gap_ns,r_gap_aux)
      }
      
      girf_inf = inf_girf_s - inf_girf_ns
      girf_ygap = ygap_girf_s - ygap_girf_ns
      girf_int = int_girf_s - int_girf_ns
      final_output_girf = rbind(c(),cbind(girf_inf,girf_ygap,girf_int))
      #colnames(final_output_girf) = c("Inf", "Ygap", "Int")
      return(final_output_girf)
      }
   
   #Test of GIRF BOOT 
   test_girfboot = GIRF_BOOT(20, mp_shock_2[2],  model_espec = list(model_inf = c("LSTR1"), model_ygap = c("LINEAR"), model_int= c("LINEAR")), state_spec = list(type = c("Indirect"), name =  c("other"), lag = 0))
   
   GIRF_SHOCK <- function(N,shock, reps){
   system.time({matrix_bootstrap = foreach(i = 1:reps, .packages = 'doParallel', .combine = cbind, .export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% {
     GIRF_BOOT(N,shock = shock,model_espec = list(model_inf = c("LSTR1"), model_ygap = c("LINEAR"), model_int= c("LINEAR")), state_spec = list(type = c("Indirect"), name =  c("other"), lag = 0))}})
   #Get Impulse Response Function For Each Variable
   GIRF_inf_t = matrix_bootstrap[,seq(1,ncol(matrix_bootstrap)-3,3)]
   GIRF_ygap_t = matrix_bootstrap[,seq(2,ncol(matrix_bootstrap)-2,3)]
   GIRF_int_t = matrix_bootstrap[,seq(3,ncol(matrix_bootstrap),3)]
   #Compute Standardized Impulse Response Function 
   SGIRF_inf_t = GIRF_inf_t/as.numeric(shock)
   SGIRF_ygap_t = GIRF_ygap_t/as.numeric(shock)
   SGIRF_int_t = GIRF_int_t/as.numeric(shock)
   
   #Compute Speed Measures (from median of standardized cumulative impulse response function):
   Speed_inf_n2_M = cumsum(rowMeans(SGIRF_inf_t))[ceiling(nrow(SGIRF_inf_t)/2)]/cumsum(rowMeans(SGIRF_inf_t))[nrow(SGIRF_inf_t)]
   Speed_inf_n4_M = cumsum(rowMeans(SGIRF_inf_t))[ceiling(nrow(SGIRF_inf_t)/4)]/cumsum(rowMeans(SGIRF_inf_t))[nrow(SGIRF_inf_t)]
   
   Speed_ygap_n2_M = cumsum(rowMeans(SGIRF_ygap_t))[ceiling(nrow(SGIRF_ygap_t)/2)]/cumsum(rowMeans(SGIRF_ygap_t))[nrow(SGIRF_ygap_t)]
   Speed_ygap_n4_M = cumsum(rowMeans(SGIRF_ygap_t))[ceiling(nrow(SGIRF_ygap_t)/4)]/cumsum(rowMeans(SGIRF_ygap_t))[nrow(SGIRF_ygap_t)]
   
   #Compute Cumulative Impulse Response and pctiles, effective, efficiency and asymmetry measures: 
   CSGIRF_inf_t = colSums(SGIRF_inf_t)
      #Measures: 
      CSGIRF_inf_t_5 = quantile(CSGIRF_inf_t, probs = 0.05)
      CSGIRF_inf_t_M = median(CSGIRF_inf_t)
      CSGIRF_inf_t_95 = quantile(CSGIRF_inf_t, probs = 0.95)
   CsGIRF_ygap_t = colSums(SGIRF_ygap_t)
   #Measures: 
      CSGIRF_ygap_t_5 = quantile(CsGIRF_ygap_t, probs = 0.05)
      CSGIRF_ygap_t_M = median(CsGIRF_ygap_t)
      CSGIRF_ygap_t_95 = quantile(CsGIRF_ygap_t, probs = 0.95)
   CsGIRF_int_t = colSums(SGIRF_int_t)
   #Measures: 
      CSGIRF_int_t_5 = quantile(CsGIRF_int_t, probs = 0.05)
      CSGIRF_int_t_M = median(CsGIRF_int_t)
      CSGIRF_int_t_95 = quantile(CsGIRF_int_t, probs = 0.95)
  #Effectiveness
      Effective_GRIF_inf = CSGIRF_inf_t_M/CSGIRF_int_t_M
      Effective_GRIF_ygap = CSGIRF_ygap_t_M/CSGIRF_int_t_M
   #Efficiency 
      Efficiency_GIRF = CSGIRF_ygap_t_M/CSGIRF_inf_t_M
   #Skewness: 
      #For inflation
      if((CSGIRF_inf_t_95-CSGIRF_inf_t_M)/(CSGIRF_inf_t_M-CSGIRF_inf_t_5) >= 1){sk_CSGIRF_inf = (CSGIRF_inf_t_95-CSGIRF_inf_t_M)/(CSGIRF_inf_t_M-CSGIRF_inf_t_5)-1} else {sk_CSGIRF_inf = 1- (CSGIRF_inf_t_M-CSGIRF_inf_t_5)/(CSGIRF_inf_t_95-CSGIRF_inf_t_M)}
      if((CSGIRF_ygap_t_95-CSGIRF_ygap_t_M)/(CSGIRF_ygap_t_M-CSGIRF_ygap_t_5) >= 1){sk_CSGIRF_ygap = (CSGIRF_ygap_t_95-CSGIRF_ygap_t_M)/(CSGIRF_ygap_t_M-CSGIRF_ygap_t_5)-1} else {sk_CSGIRF_ygap = 1-(CSGIRF_ygap_t_M-CSGIRF_ygap_t_5)/(CSGIRF_ygap_t_95-CSGIRF_ygap_t_M)}
      if((CSGIRF_int_t_95-CSGIRF_int_t_M)/(CSGIRF_int_t_M-CSGIRF_int_t_5) >= 1){sk_CSGIRF_int = (CSGIRF_int_t_95-CSGIRF_int_t_M)/(CSGIRF_int_t_M-CSGIRF_int_t_5)-1} else {sk_CSGIRF_int = 1- (CSGIRF_int_t_M-CSGIRF_int_t_5)/(CSGIRF_int_t_95-CSGIRF_int_t_M)}
   return(rbind(c(), cbind(CSGIRF_inf_t_5,CSGIRF_inf_t_M,CSGIRF_inf_t_95,CSGIRF_ygap_t_5,CSGIRF_ygap_t_M,CSGIRF_ygap_t_95,CSGIRF_int_t_5,CSGIRF_int_t_M,CSGIRF_int_t_95,Effective_GRIF_inf,Effective_GRIF_ygap,Efficiency_GIRF,sk_CSGIRF_inf,sk_CSGIRF_ygap,sk_CSGIRF_int, Speed_inf_n2_M,Speed_inf_n4_M, Speed_ygap_n2_M, Speed_ygap_n4_M)))
   } 
   
   #Test GIRF SHOCK
   #test_girfshock = GIRF_SHOCK(20,mp_shock_2[5], reps = 1000)
   
   system.time({hola = foreach(i = 1:2, .packages = c('doParallel'), .combine = cbind, .export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% GIRF_SHOCK(N = 20, mp_shock_2[i], reps = 1000)})

   
   #test_girft = GIRF_T(1)
   
   #Data for propagation:    
   data_girf_inf = data_str_inf_model[-1, -2]
   data_girf_outputgap = data_str_outputgap_model[-1,-2]
   data_girf_interest = data_str_interest_model[-1,-2]
   real_neutral_rate = data_xts$real_natural_rate_col["2003-03-01/2019-12-01"]
   
   # Start the clock!
   ptm <- proc.time()
   system.time({GIRF_FINAL = foreach(x = 1:7, .combine = rbind, .packages = c('doParallel', 'minpack.lm'),.export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% 
      {
         data_girf_inf = data_girf_inf[x:nrow(data_girf_inf),]
         data_girf_outputgap = data_girf_outputgap[x:nrow(data_girf_outputgap),]
         data_girf_interest = data_girf_interest[x:nrow(data_girf_interest),]
         real_neutral_rate = real_neutral_rate[x:nrow(real_neutral_rate),]
         state_var_inf_girf = state_var_inf_girf[x:length(state_var_inf_girf)]
         #state_var_ygap_girf = state_var_ygap_girf[x:length(state_var_ygap_girf)]
         #state_var_int_girf = state_var_int_girf[x:length(state_var_int_girf)] 
         foreach(i = 1:101, .packages = c('doParallel'), .combine = cbind, .export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% {GIRF_SHOCK(20, mp_shock_2[i], reps = 1000)}         
      
   }})    

   system.time({GIRF_FINAL_2 = foreach(x = 8:14, .combine = rbind, .packages = c('doParallel', 'minpack.lm'),.export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% 
      {
         data_girf_inf = data_girf_inf[x:nrow(data_girf_inf),]
         data_girf_outputgap = data_girf_outputgap[x:nrow(data_girf_outputgap),]
         data_girf_interest = data_girf_interest[x:nrow(data_girf_interest),]
         real_neutral_rate = real_neutral_rate[x:nrow(real_neutral_rate),]
         state_var_inf_girf = state_var_inf_girf[x:length(state_var_inf_girf)]
         #state_var_ygap_girf = state_var_ygap_girf[x:length(state_var_ygap_girf)]
         #state_var_int_girf = state_var_int_girf[x:length(state_var_int_girf)] 
         foreach(i = 1:101, .packages = c('doParallel'), .combine = cbind, .export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% {GIRF_SHOCK(20,mp_shock_2[i], reps = 1000)}         
         
      }})  
 
  
   system.time({GIRF_FINAL_3 = foreach(x = 15:21, .combine = rbind, .packages = c('doParallel', 'minpack.lm'),.export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% 
      {
         data_girf_inf = data_girf_inf[x:nrow(data_girf_inf),]
         data_girf_outputgap = data_girf_outputgap[x:nrow(data_girf_outputgap),]
         data_girf_interest = data_girf_interest[x:nrow(data_girf_interest),]
         real_neutral_rate = real_neutral_rate[x:nrow(real_neutral_rate),]
         state_var_inf_girf = state_var_inf_girf[x:length(state_var_inf_girf)]
         #state_var_ygap_girf = state_var_ygap_girf[x:length(state_var_ygap_girf)]
         #state_var_int_girf = state_var_int_girf[x:length(state_var_int_girf)] 
         foreach(i = 1:101, .packages = c('doParallel'), .combine = cbind, .export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% {GIRF_SHOCK(20,mp_shock_2[i], reps = 1000)}         
         
      }})  
   
   
   system.time({GIRF_FINAL_4 = foreach(x = 22:28, .combine = rbind, .packages = c('doParallel', 'minpack.lm'),.export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% 
      {
         data_girf_inf = data_girf_inf[x:nrow(data_girf_inf),]
         data_girf_outputgap = data_girf_outputgap[x:nrow(data_girf_outputgap),]
         data_girf_interest = data_girf_interest[x:nrow(data_girf_interest),]
         real_neutral_rate = real_neutral_rate[x:nrow(real_neutral_rate),]
         state_var_inf_girf = state_var_inf_girf[x:length(state_var_inf_girf)]
         #state_var_ygap_girf = state_var_ygap_girf[x:length(state_var_ygap_girf)]
         #state_var_int_girf = state_var_int_girf[x:length(state_var_int_girf)] 
         foreach(i = 1:101, .packages = c('doParallel'), .combine = cbind, .export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% {GIRF_SHOCK(20,mp_shock_2[i], reps = 1000)}         
         
      }})
      
   system.time({GIRF_FINAL_5 = foreach(x = 29:35, .combine = rbind, .packages = c('doParallel', 'minpack.lm'),.export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% 
      {
         data_girf_inf = data_girf_inf[x:nrow(data_girf_inf),]
         data_girf_outputgap = data_girf_outputgap[x:nrow(data_girf_outputgap),]
         data_girf_interest = data_girf_interest[x:nrow(data_girf_interest),]
         real_neutral_rate = real_neutral_rate[x:nrow(real_neutral_rate),]
         state_var_inf_girf = state_var_inf_girf[x:length(state_var_inf_girf)]
         #state_var_ygap_girf = state_var_ygap_girf[x:length(state_var_ygap_girf)]
         #state_var_int_girf = state_var_int_girf[x:length(state_var_int_girf)] 
         foreach(i = 1:101, .packages = c('doParallel'), .combine = cbind, .export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% {GIRF_SHOCK(20,mp_shock_2[i], reps = 1000)}         
         
      }})
      
      system.time({GIRF_FINAL_6 = foreach(x = 36:42, .combine = rbind, .packages = c('doParallel', 'minpack.lm'),.export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% 
         {
            data_girf_inf = data_girf_inf[x:nrow(data_girf_inf),]
            data_girf_outputgap = data_girf_outputgap[x:nrow(data_girf_outputgap),]
            data_girf_interest = data_girf_interest[x:nrow(data_girf_interest),]
            real_neutral_rate = real_neutral_rate[x:nrow(real_neutral_rate),]
            state_var_inf_girf = state_var_inf_girf[x:length(state_var_inf_girf)]
            #state_var_ygap_girf = state_var_ygap_girf[x:length(state_var_ygap_girf)]
            #state_var_int_girf = state_var_int_girf[x:length(state_var_int_girf)] 
            foreach(i = 1:101, .packages = c('doParallel'), .combine = cbind, .export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% {GIRF_SHOCK(20,mp_shock_2[i], reps = 1000)}         
            
         }})
   
   system.time({GIRF_FINAL_7 = foreach(x = 43:47, .combine = rbind, .packages = c('doParallel', 'minpack.lm'),.export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% 
      {
         data_girf_inf = data_girf_inf[x:nrow(data_girf_inf),]
         data_girf_outputgap = data_girf_outputgap[x:nrow(data_girf_outputgap),]
         data_girf_interest = data_girf_interest[x:nrow(data_girf_interest),]
         real_neutral_rate = real_neutral_rate[x:nrow(real_neutral_rate),]
         state_var_inf_girf = state_var_inf_girf[x:length(state_var_inf_girf)]
         #state_var_ygap_girf = state_var_ygap_girf[x:length(state_var_ygap_girf)]
         #state_var_int_girf = state_var_int_girf[x:length(state_var_int_girf)] 
         foreach(i = 1:101, .packages = c('doParallel'), .combine = cbind, .export = c("state_var_inf_girf", "state_var_int_girf", "state_var_ygap_girf")) %dopar% {GIRF_SHOCK(20,mp_shock_2[i], reps = 1000)}         
         
      }})
   # Stop the clock
   proc.time() - ptm
   
   
   #GIRF_state_iedp = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_services = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_vivicr = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_invgdp = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_brent = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_tasaext = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_vencida = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_iicb = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_inf_ext = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_deuda_ext = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_carterasecun = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_comercialicm = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_bonospasiv = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_hurtog = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_offenfarc = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_IHHcdt = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_IHHcomer = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_IHHcartera = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_tasacomerdes = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_totalicm = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_cubrtotal = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_consumcubr = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_agricul = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_PIBd = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   GIRF_state_IICV = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_ICE = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_outputgap_l1 = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   #GIRF_state_outputgap = rbind(GIRF_FINAL, GIRF_FINAL_2, GIRF_FINAL_3,GIRF_FINAL_4,GIRF_FINAL_5,GIRF_FINAL_6,GIRF_FINAL_7)
   
   #Data set for State and Sign Asymmetries Regressions and Tests: 
      #Sign and Size Asymmetry measures
      length(mp_shock_2)
         #Sign :
         positive_shock_sign = mp_shock_2
         positive_shock_sign[which(positive_shock_sign >= 0)] <- 1
         positive_shock_sign[which(positive_shock_sign < 0)] <- 0
         negative_shock_sign = mp_shock_2
         negative_shock_sign[which(negative_shock_sign >= 0)] <- 0
         negative_shock_sign[which(negative_shock_sign < 0)] <- 1
         
         #Size: 
         shock_size_positive = mp_shock_2
         shock_size_positive[which(shock_size_positive < 0)] <- 0
         shock_size_negative = mp_shock_2
         shock_size_negative[which(shock_size_negative >= 0)] <- 0
   
   #Measures of interest:       
   data_state_sign = colMeans(GIRF_state_IICV) #Change data set
   
   CSGIRF_inf_M = data_state_sign[seq(2,length(data_state_sign),19)]
   Effective_inf = data_state_sign[seq(10,length(data_state_sign),19)]
   Skew_inf = data_state_sign[seq(13,length(data_state_sign),19)]
   Speed_2_inf = data_state_sign[seq(16,length(data_state_sign),19)]
   Speed_4_inf = data_state_sign[seq(17,length(data_state_sign),19)]
   
   CSGIRF_ygap_M = data_state_sign[seq(5,length(data_state_sign),19)]
   Effective_ygap = data_state_sign[seq(11,length(data_state_sign),19)]
   Skew_ygap = data_state_sign[seq(14,length(data_state_sign),19)]
   Speed_2_ygap = data_state_sign[seq(18,length(data_state_sign),19)]
   Speed_4_ygap = data_state_sign[seq(19,length(data_state_sign),19)]
   
   CSGIRF_interest_M = data_state_sign[seq(8,length(data_state_sign),19)]
   Skew_interest = data_state_sign[seq(15,length(data_state_sign),19)]
   
   Efficiency_GIRF = data_state_sign[seq(12,length(data_state_sign),19)]
   
   #Sign and Size Asymmetry regressions (OLS)
   
   #Inflation Equation Size and Sign Asymmetry Tests 
   lm_inf_ssa = lm(CSGIRF_inf_M ~ positive_shock_sign + negative_shock_sign + shock_size_positive+shock_size_negative-1)
   summary(lm_inf_ssa)   
      #Sign Asymmetry OVerall
      wald_1 = wald.test(vcov(lm_inf_ssa), coef(lm_inf_ssa), c(1,2,3,4), H0 = c(coef(lm_inf_ssa)[2],coef(lm_inf_ssa)[1],coef(lm_inf_ssa)[4],coef(lm_inf_ssa)[3]),df = lm_inf_ssa$df.residual)
      #Sign Asymmetry
      wald_2 = wald.test(vcov(lm_inf_ssa), coef(lm_inf_ssa), c(1,2), H0 = c(coef(lm_inf_ssa)[2],coef(lm_inf_ssa)[1]), df = lm_inf_ssa$df.residual)
      #Size Asymmetry, if it exists, is equal for positive and negative shocks 
      wald_3 = wald.test(vcov(lm_inf_ssa), coef(lm_inf_ssa), c(3,4), H0 = c(coef(lm_inf_ssa)[4],coef(lm_inf_ssa)[3]), df = lm_inf_ssa$df.residual)
      #Size Asymmetry
      wald_4 = wald.test(vcov(lm_inf_ssa), coef(lm_inf_ssa), c(3,4), H0 = c(0,0), df = lm_inf_ssa$df.residual)

   salida_wald_ssa  = cbind(wald_1$result$Ftest[4],wald_1$result$chi2[3],wald_2$result$Ftest[4],wald_2$result$chi2[3],wald_3$result$Ftest[4],wald_3$result$chi2[3],wald_4$result$Ftest[4],wald_4$result$chi2[3])
   salida_wald_ssa   
   
   lm_infef_ssa = lm(Effective_inf ~ positive_shock_sign + negative_shock_sign + shock_size_positive+shock_size_negative-1)
   summary(lm_infef_ssa)   
      #Sign Asymmetry OVerall
      wald_1 = wald.test(vcov(lm_infef_ssa), coef(lm_infef_ssa), c(1,2,3,4), H0 = c(coef(lm_infef_ssa)[2],coef(lm_infef_ssa)[1],coef(lm_infef_ssa)[4],coef(lm_infef_ssa)[3]), df = lm_infef_ssa$df.residual)
      #Sign Asymetry
      wald_2 = wald.test(vcov(lm_infef_ssa), coef(lm_infef_ssa), c(1,2), H0 = c(coef(lm_infef_ssa)[2],coef(lm_infef_ssa)[1]), df = lm_infef_ssa$df.residual)
      #Size Asymmetry, if it exists, is equal for positive and negative shocks 
      wald_3 = wald.test(vcov(lm_infef_ssa), coef(lm_infef_ssa), c(3,4), H0 = c(coef(lm_infef_ssa)[4],coef(lm_infef_ssa)[3]), df = lm_infef_ssa$df.residual)
      #Size Asymmetry
      wald_4 = wald.test(vcov(lm_infef_ssa), coef(lm_infef_ssa), c(3,4), H0 = c(0,0), df = lm_infef_ssa$df.residual)
      
   salida_wald_ssa  = cbind(wald_1$result$Ftest[4],wald_1$result$chi2[3],wald_2$result$Ftest[4],wald_2$result$chi2[3],wald_3$result$Ftest[4],wald_3$result$chi2[3],wald_4$result$Ftest[4],wald_4$result$chi2[3])
   salida_wald_ssa   
            
   lm_infskew_ssa = lm(Skew_inf ~ positive_shock_sign + negative_shock_sign + shock_size_positive+shock_size_negative-1)
   summary(lm_infskew_ssa)   
      #Non Zero Skew
      wald_1 = wald.test(vcov(lm_infskew_ssa), coef(lm_infskew_ssa), c(1,2), H0 = c(0,0), df = lm_infskew_ssa$df.residual)
      #Skew Sign Asymetry
      wald_2 = wald.test(vcov(lm_infskew_ssa), coef(lm_infskew_ssa), c(1,2), H0 = c(coef(lm_infskew_ssa)[2],coef(lm_infskew_ssa)[1]), df = lm_infskew_ssa$df.residual)
      wald_3 = wald.test(vcov(lm_infskew_ssa), coef(lm_infskew_ssa), c(1,2), H0 = c(-coef(lm_infskew_ssa)[2],-coef(lm_infskew_ssa)[1]), df = lm_infskew_ssa$df.residual)
      #Size Asymmetry
      wald_4 = wald.test(vcov(lm_infskew_ssa), coef(lm_infskew_ssa), c(3,4), H0 = c(0,0), df = lm_infskew_ssa$df.residual)
      #Size and Sign Asymetry
      wald_5 = wald.test(vcov(lm_infskew_ssa), coef(lm_infskew_ssa), c(3,4), H0 = c(coef(lm_infskew_ssa)[4],coef(lm_infskew_ssa)[3]), df = lm_infskew_ssa$df.residual)
      wald_6 = wald.test(vcov(lm_infskew_ssa), coef(lm_infskew_ssa), c(3,4), H0 = c(-coef(lm_infskew_ssa)[4],-coef(lm_infskew_ssa)[3]), df = lm_infskew_ssa$df.residual)
      
      salida_wald_ssa  = cbind(wald_1$result$Ftest[4],wald_1$result$chi2[3],wald_2$result$Ftest[4],wald_2$result$chi2[3],wald_3$result$Ftest[4],wald_3$result$chi2[3],wald_4$result$Ftest[4],wald_4$result$chi2[3],wald_5$result$Ftest[4],wald_5$result$chi2[3],wald_6$result$Ftest[4],wald_6$result$chi2[3])
      salida_wald_ssa
               
   lm_infspeed4_ssa = lm(Speed_4_inf ~ positive_shock_sign + negative_shock_sign + shock_size_positive+shock_size_negative-1)
   summary(lm_infspeed4_ssa)   
      #Sign Asymmetry OVerall
      wald_1 = wald.test(vcov(lm_infspeed4_ssa), coef(lm_infspeed4_ssa), c(1,2,3,4), H0 = c(coef(lm_infspeed4_ssa)[2],coef(lm_infspeed4_ssa)[1],coef(lm_infspeed4_ssa)[4],coef(lm_infspeed4_ssa)[3]), df = lm_infspeed4_ssa$df.residual)
      #Sign Asymmetry
      wald_2 = wald.test(vcov(lm_infspeed4_ssa), coef(lm_infspeed4_ssa), c(1,2), H0 = c(coef(lm_infspeed4_ssa)[2],coef(lm_infspeed4_ssa)[1]), df = lm_infspeed4_ssa$df.residual)
      #Size Asymmetry, if it exists, is equal for positive and negative shocks 
      wald_3 = wald.test(vcov(lm_infspeed4_ssa), coef(lm_infspeed4_ssa), c(3,4), H0 = c(coef(lm_infspeed4_ssa)[4],coef(lm_infspeed4_ssa)[3]), df = lm_infspeed4_ssa$df.residual)
      #Size Asymmetry
      wald_4 = wald.test(vcov(lm_infspeed4_ssa), coef(lm_infspeed4_ssa), c(3,4), H0 = c(0,0), df = lm_infspeed4_ssa$df.residual)
   
   salida_wald_ssa  = cbind(wald_1$result$Ftest[4],wald_1$result$chi2[3],wald_2$result$Ftest[4],wald_2$result$chi2[3],wald_3$result$Ftest[4],wald_3$result$chi2[3],wald_4$result$Ftest[4],wald_4$result$chi2[3])
   salida_wald_ssa
   
   #Output Gap Sign and Size Asymmetry Test
         
   lm_ygap_ssa = lm(CSGIRF_ygap_M ~ positive_shock_sign + negative_shock_sign + shock_size_positive+shock_size_negative-1)   
   summary(lm_ygap_ssa)   
      #Sign Asymmetry OVerall
      wald_1 =  wald.test(vcov(lm_ygap_ssa), coef(lm_ygap_ssa), c(1,2,3,4), H0 = c(coef(lm_ygap_ssa)[2],coef(lm_ygap_ssa)[1],coef(lm_ygap_ssa)[4],coef(lm_ygap_ssa)[3]), df = lm_ygap_ssa$df.residual)
      #Sign Asymmetry   
      wald_2 =  wald.test(vcov(lm_ygap_ssa), coef(lm_ygap_ssa), c(1,2), H0 = c(coef(lm_ygap_ssa)[2],coef(lm_ygap_ssa)[1]), df = lm_ygap_ssa$df.residual)
      #Size Asymmetry, if it exists, is equal for positive and negative shocks 
      wald_3 = wald.test(vcov(lm_ygap_ssa), coef(lm_ygap_ssa), c(3,4), H0 = c(coef(lm_ygap_ssa)[4],coef(lm_ygap_ssa)[3]), df = lm_ygap_ssa$df.residual)
      #Size Asymmetry
      wald_4 = wald.test(vcov(lm_ygap_ssa), coef(lm_ygap_ssa), c(3,4), H0 = c(0,0), df = lm_ygap_ssa$df.residual)

   salida_wald_ssa  = cbind(wald_1$result$Ftest[4],wald_1$result$chi2[3],wald_2$result$Ftest[4],wald_2$result$chi2[3],wald_3$result$Ftest[4],wald_3$result$chi2[3],wald_4$result$Ftest[4],wald_4$result$chi2[3])
   salida_wald_ssa
      
   lm_ygapef_ssa = lm(Effective_ygap ~ positive_shock_sign + negative_shock_sign + shock_size_positive+shock_size_negative-1)   
   summary(lm_ygapef_ssa)
      #Sign Asymmetry OVerall
      wald_1 = wald.test(vcov(lm_ygapef_ssa), coef(lm_ygapef_ssa), c(1,2,3,4), H0 = c(coef(lm_ygapef_ssa)[2],coef(lm_ygapef_ssa)[1],coef(lm_ygapef_ssa)[4],coef(lm_ygapef_ssa)[3]), df = lm_ygapef_ssa$df.residual)
      #Sign Asymmetry   
      wald_2 = wald.test(vcov(lm_ygapef_ssa), coef(lm_ygapef_ssa), c(1,2), H0 = c(coef(lm_ygapef_ssa)[2],coef(lm_ygapef_ssa)[1]), df = lm_ygapef_ssa$df.residual)
      #Size Asymmetry, if it exists, is equal for positive and negative shocks 
      wald_3 = wald.test(vcov(lm_ygapef_ssa), coef(lm_ygapef_ssa), c(3,4), H0 = c(coef(lm_ygapef_ssa)[4],coef(lm_ygapef_ssa)[3]), df = lm_ygapef_ssa$df.residual)
      #Size Asymmetry
      wald_4 = wald.test(vcov(lm_ygapef_ssa), coef(lm_ygapef_ssa), c(3,4), H0 = c(0,0), df = lm_ygapef_ssa$df.residual)
      
      salida_wald_ssa  = cbind(wald_1$result$Ftest[4],wald_1$result$chi2[3],wald_2$result$Ftest[4],wald_2$result$chi2[3],wald_3$result$Ftest[4],wald_3$result$chi2[3],wald_4$result$Ftest[4],wald_4$result$chi2[3])
      salida_wald_ssa
      
   lm_ygapskew_ssa = lm(Skew_ygap ~ positive_shock_sign + negative_shock_sign + shock_size_positive+shock_size_negative-1)
   summary(lm_ygapskew_ssa)   
      #Non Zero Skew
      wald_1 = wald.test(vcov(lm_ygapskew_ssa), coef(lm_ygapskew_ssa), c(1,2), H0 = c(0,0), df = lm_ygapskew_ssa$df.residual)
      #Skew Sign Asymetry
      wald_2 =wald.test(vcov(lm_ygapskew_ssa), coef(lm_ygapskew_ssa), c(1,2), H0 = c(coef(lm_ygapskew_ssa)[2],coef(lm_ygapskew_ssa)[1]), df = lm_ygapskew_ssa$df.residual)
      wald_3 =wald.test(vcov(lm_ygapskew_ssa), coef(lm_ygapskew_ssa), c(1,2), H0 = c(-coef(lm_ygapskew_ssa)[2],-coef(lm_ygapskew_ssa)[1]), df = lm_ygapskew_ssa$df.residual)
      #Size Asymmetry
      wald_4 =wald.test(vcov(lm_ygapskew_ssa), coef(lm_ygapskew_ssa), c(3,4), H0 = c(0,0), df = lm_infskew_ssa$df.residual)
      #Size and Sign Asymetry
      wald_5 =   wald.test(vcov(lm_ygapskew_ssa), coef(lm_ygapskew_ssa), c(3,4), H0 = c(coef(lm_ygapskew_ssa)[4],coef(lm_ygapskew_ssa)[3]), df = lm_ygapskew_ssa$df.residual)
      wald_6 =   wald.test(vcov(lm_ygapskew_ssa), coef(lm_ygapskew_ssa), c(3,4), H0 = c(-coef(lm_ygapskew_ssa)[4],-coef(lm_ygapskew_ssa)[3]), df = lm_ygapskew_ssa$df.residual)
   
      salida_wald_ssa  = cbind(wald_1$result$Ftest[4],wald_1$result$chi2[3],wald_2$result$Ftest[4],wald_2$result$chi2[3],wald_3$result$Ftest[4],wald_3$result$chi2[3],wald_4$result$Ftest[4],wald_4$result$chi2[3],wald_5$result$Ftest[4],wald_5$result$chi2[3],wald_6$result$Ftest[4],wald_6$result$chi2[3])
      salida_wald_ssa
      
   lm_ygapspeed4_ssa = lm(Speed_4_ygap ~ positive_shock_sign + negative_shock_sign + shock_size_positive+shock_size_negative-1)   
   summary(lm_ygapspeed4_ssa)
      #Sign Asymmetry OVerall
      wald_1 = wald.test(vcov(lm_ygapspeed4_ssa), coef(lm_ygapspeed4_ssa), c(1,2,3,4), H0 = c(coef(lm_ygapspeed4_ssa)[2],coef(lm_ygapspeed4_ssa)[1],coef(lm_ygapspeed4_ssa)[4],coef(lm_ygapspeed4_ssa)[3]), df = lm_ygapspeed4_ssa$df.residual)
      #Sign Asymmetry   
      wald_2 = wald.test(vcov(lm_ygapspeed4_ssa), coef(lm_ygapspeed4_ssa), c(1,2), H0 = c(coef(lm_ygapspeed4_ssa)[2],coef(lm_ygapspeed4_ssa)[1]), df = lm_ygapspeed4_ssa$df.residual)
      #Size Asymmetry, if it exists, is equal for positive and negative shocks 
      wald_3 = wald.test(vcov(lm_ygapspeed4_ssa), coef(lm_ygapspeed4_ssa), c(3,4), H0 = c(coef(lm_ygapspeed4_ssa)[4],coef(lm_ygapspeed4_ssa)[3]), df = lm_ygapspeed4_ssa$df.residual)
      #Size Asymmetry
      wald_4 = wald.test(vcov(lm_ygapspeed4_ssa), coef(lm_ygapspeed4_ssa), c(3,4), H0 = c(0,0), df = lm_ygapspeed4_ssa$df.residual)
      
      salida_wald_ssa  = cbind(wald_1$result$Ftest[4],wald_1$result$chi2[3],wald_2$result$Ftest[4],wald_2$result$chi2[3],wald_3$result$Ftest[4],wald_3$result$chi2[3],wald_4$result$Ftest[4],wald_4$result$chi2[3])
      salida_wald_ssa    
      
   #Interest Rate Equtaion Sign and Size Asymmetry 
           
   lm_interest_ssa = lm(CSGIRF_interest_M ~ positive_shock_sign + negative_shock_sign + shock_size_positive+shock_size_negative-1)   
    summary(lm_interest_ssa) 
      #Sign Asymmetry OVerall
      wald_1 = wald.test(vcov(lm_interest_ssa), coef(lm_interest_ssa), c(1,2,3,4), H0 = c(coef(lm_interest_ssa)[2],coef(lm_interest_ssa)[1],coef(lm_interest_ssa)[4],coef(lm_interest_ssa)[3]), df = lm_interest_ssa$df.residual)
      #Sign Asymmetry
      wald_2 = wald.test(vcov(lm_interest_ssa), coef(lm_interest_ssa), c(1,2), H0 = c(coef(lm_interest_ssa)[2],coef(lm_interest_ssa)[1]), df = lm_interest_ssa$df.residual)
      #Size Asymmetry, if it exists, is equal for positive and negative shocks 
      wald_3 = wald.test(vcov(lm_interest_ssa), coef(lm_interest_ssa), c(3,4), H0 = c(coef(lm_interest_ssa)[4],coef(lm_interest_ssa)[3]), df = lm_interest_ssa$df.residual)
      #Size Asymmetry
      wald_4 = wald.test(vcov(lm_interest_ssa), coef(lm_interest_ssa), c(3,4), H0 = c(0,0), df = lm_interest_ssa$df.residual)   
      
      salida_wald_ssa  = cbind(wald_1$result$Ftest[4],wald_1$result$chi2[3],wald_2$result$Ftest[4],wald_2$result$chi2[3],wald_3$result$Ftest[4],wald_3$result$chi2[3],wald_4$result$Ftest[4],wald_4$result$chi2[3])
      salida_wald_ssa
      
   lm_interestskew_ssa = lm(Skew_interest ~ positive_shock_sign + negative_shock_sign + shock_size_positive+shock_size_negative-1)
   summary(lm_interestskew_ssa)   
      #Non Zero Skew
      wald_1 = wald.test(vcov(lm_interestskew_ssa), coef(lm_interestskew_ssa), c(1,2), H0 = c(0,0), df = lm_interestskew_ssa$df.residual)
      #Skew Sign Asymetry
      wald_2 = wald.test(vcov(lm_interestskew_ssa), coef(lm_interestskew_ssa), c(1,2), H0 = c(coef(lm_interestskew_ssa)[2],coef(lm_interestskew_ssa)[1]), df = lm_interestskew_ssa$df.residual)
      wald_3 = wald.test(vcov(lm_interestskew_ssa), coef(lm_interestskew_ssa), c(1,2), H0 = c(-coef(lm_interestskew_ssa)[2],-coef(lm_interestskew_ssa)[1]), df = lm_interestskew_ssa$df.residual)
      #Size Asymmetry
      wald_4 = wald.test(vcov(lm_interestskew_ssa), coef(lm_interestskew_ssa), c(3,4), H0 = c(0,0), df = lm_infskew_ssa$df.residual)
      #Size and Sign Asymetry
      wald_5 = wald.test(vcov(lm_interestskew_ssa), coef(lm_interestskew_ssa), c(3,4), H0 = c(coef(lm_interestskew_ssa)[4],coef(lm_interestskew_ssa)[3]), df = lm_interestskew_ssa$df.residual)
      wald_6 = wald.test(vcov(lm_interestskew_ssa), coef(lm_interestskew_ssa), c(3,4), H0 = c(-coef(lm_interestskew_ssa)[4],-coef(lm_interestskew_ssa)[3]), df = lm_interestskew_ssa$df.residual)
      
      salida_wald_ssa  = cbind(wald_1$result$Ftest[4],wald_1$result$chi2[3],wald_2$result$Ftest[4],wald_2$result$chi2[3],wald_3$result$Ftest[4],wald_3$result$chi2[3],wald_4$result$Ftest[4],wald_4$result$chi2[3],wald_5$result$Ftest[4],wald_5$result$chi2[3],wald_6$result$Ftest[4],wald_6$result$chi2[3])
      salida_wald_ssa
      
   #Efficiency sign and size asymmetry 
   lm_efficiency_ssa = lm(Efficiency_GIRF ~ positive_shock_sign + negative_shock_sign + shock_size_positive+shock_size_negative-1)
   summary(lm_efficiency_ssa)   
      #Sign Asymmetry OVerall
      wald_1 = wald.test(vcov(lm_efficiency_ssa), coef(lm_efficiency_ssa), c(1,2,3,4), H0 = c(coef(lm_efficiency_ssa)[2],coef(lm_efficiency_ssa)[1],coef(lm_efficiency_ssa)[4],coef(lm_efficiency_ssa)[3]), df = lm_efficiency_ssa$df.residual)
      #Sign Asymmetry
      wald_2 = wald.test(vcov(lm_efficiency_ssa), coef(lm_efficiency_ssa), c(1,2), H0 = c(coef(lm_efficiency_ssa)[2],coef(lm_efficiency_ssa)[1]), df = lm_efficiency_ssa$df.residual)
      #Size Asymmetry, if it exists, is equal for positive and negative shocks 
      wald_3 = wald.test(vcov(lm_efficiency_ssa), coef(lm_efficiency_ssa), c(3,4), H0 = c(coef(lm_efficiency_ssa)[4],coef(lm_efficiency_ssa)[3]), df = lm_efficiency_ssa$df.residual)
      #Size Asymmetry
      wald_4 = wald.test(vcov(lm_efficiency_ssa), coef(lm_efficiency_ssa), c(3,4), H0 = c(0,0), df = lm_efficiency_ssa$df.residual)   
   
      salida_wald_ssa  = cbind(wald_1$result$Ftest[4],wald_1$result$chi2[3],wald_2$result$Ftest[4],wald_2$result$chi2[3],wald_3$result$Ftest[4],wald_3$result$chi2[3],wald_4$result$Ftest[4],wald_4$result$chi2[3])
      salida_wald_ssa
      
   salida_regresiones_inf = cbind(summary(lm_inf_ssa)$coefficients[,1],summary(lm_inf_ssa)$coefficients[,4],
                              summary(lm_infef_ssa)$coefficients[,1],summary(lm_infef_ssa)$coefficients[,4],
                              summary(lm_infskew_ssa)$coefficients[,1],summary(lm_infskew_ssa)$coefficients[,4],
                              summary(lm_infspeed4_ssa)$coefficients[,1],summary(lm_infspeed4_ssa)$coefficients[,4]) 
   salida_regresiones_ygap = cbind(summary(lm_ygap_ssa)$coefficients[,1],summary(lm_ygap_ssa)$coefficients[,4],
                                   summary(lm_ygapef_ssa)$coefficients[,1],summary(lm_ygapef_ssa)$coefficients[,4],
                                   summary(lm_ygapskew_ssa)$coefficients[,1],summary(lm_ygapskew_ssa)$coefficients[,4],
                                   summary(lm_ygapspeed4_ssa)$coefficients[,1],summary(lm_ygapspeed4_ssa)$coefficients[,4])

   
         
   #Data for State Asymmetry (change here state variable)
   state_variable_lm = NULL 
   for(i in seq(1,nrow(LSTR1_INF$IICV$Transition_function)-21,1)){ #The change of state variable should go here
     aux = mean(LSTR1_INF$IICV$Transition_function[i:eval(i+20),3]) 
     state_variable_lm[i] = aux
   }
   plot(state_variable_lm, type = 'l')
   
   GIRF_data_aux = GIRF_state_outputgap #Change data set
      #Import GIRF and related measures of 10%, 25%, 50%, 75% and 90% 
      Import_GIRF_measures = cbind(GIRF_data_aux[,eval(10*19+1):eval(11*19)],GIRF_data_aux[,eval(25*19+1):eval(26*19)],GIRF_data_aux[,eval(50*19+1):eval(51*19)],GIRF_data_aux[,eval(75*19+1):eval(76*19)],GIRF_data_aux[,eval(90*19+1):eval(91*19)])
      #write.csv(Import_GIRF_measures, "GIRF_measures_outputgapl1.csv")
   
      GIRF_M_import = Import_GIRF_measures[,c(2,21,40,59,78,5,24,43,62,81)]
      GIRF_speed_import = Import_GIRF_measures[,c(17,36,55,74,93,19,38,57,76,95)]
      GIRF_efficiency_import = Import_GIRF_measures[,c(12,31,50,69,88)]
      GIRF_effective_import = Import_GIRF_measures[,c(10,29,48,67,86,11,30,49,68,87)]
      GIRF_risk_import = Import_GIRF_measures[,c(13,32,51,70,89,14,33,52,71,90)]
      
      write.csv(GIRF_M_import, 'GIRF_M_inf_ext.csv')
      write.csv(GIRF_speed_import, 'GIRF_speed_inf_ext.csv')
      write.csv(GIRF_efficiency_import, 'GIRF_efficiency_outputgap.csv')
      write.csv(GIRF_effective_import, 'GIRF_effective_outputgap.csv')
      write.csv(GIRF_risk_import, 'GIRF_risk_inf_ext.csv')
 
# ---------------------------- Regresiones Rolling -------------------------- #      
      #Choque positivo
      GIRF_data_aux_ps_90 = GIRF_data_aux[,eval(90*19+1):eval(91*19)]
      colnames_salida_state = colnames(GIRF_data_aux_ps_90)
      salida_state_ps_90_rolling = NULL
      for(j in 0:eval(nrow(GIRF_data_aux_ps_90)-30)){
         salida_state_ps_90 = NULL
         for(i in 1:ncol(GIRF_data_aux_ps_90)){
            model_aux = summary(lm(GIRF_data_aux_ps_90[eval(1+j):eval(30+j),i]~ 1+state_variable_lm[eval(1+j):eval(30+j)]))
            model_aux
            beta = model_aux$coefficients[2,1]
            std = model_aux$coefficients[2,2]
            pval = model_aux$coefficients[2,4]
            salida_aux = cbind(beta,std, pval)
            salida_state_ps_90 = rbind(salida_state_ps_90, salida_aux)
         }
         rownames(salida_state_ps_90) = colnames_salida_state
         aux = salida_state_ps_90
         salida_state_ps_90_rolling = cbind(salida_state_ps_90_rolling,aux)
      }
      salida_state_ps_90_rolling[c(2,5,10,11,12,13,14,17,19),]
      write.csv(salida_state_ps_90_rolling[c(2,5,10,11,12,13,14,17,19),], 'betas_ict_ps_90_rolling.csv')
      
      # Choque negativo 
      GIRF_data_aux_ns_10 = GIRF_data_aux[,eval(10*19+1):eval(11*19)]
      colnames_salida_state = colnames(GIRF_data_aux_ns_10)
      salida_state_ns_10_rolling = NULL
      for(j in 0:eval(nrow(GIRF_data_aux_ns_10)-30)){
         salida_state_ns_10 = NULL
         for(i in 1:ncol(GIRF_data_aux_ns_10)){
            model_aux = summary(lm(GIRF_data_aux_ns_10[eval(1+j):eval(30+j),i]~ 1+state_variable_lm[eval(1+j):eval(30+j)]))
            model_aux
            beta = model_aux$coefficients[2,1]
            std = model_aux$coefficients[2,2]
            pval = model_aux$coefficients[2,4]
            salida_aux = cbind(beta, std, pval)
            salida_state_ns_10 = rbind(salida_state_ns_10, salida_aux)
         }
         rownames(salida_state_ns_10) = colnames_salida_state
         aux = salida_state_ns_10
         salida_state_ns_10_rolling = cbind(salida_state_ns_10_rolling, aux)
      }
      write.csv(salida_state_ns_10_rolling[c(2,5,10,11,12,13,14,17,19),], 'betas_ict_ps_10_rolling.csv')
      
      
# ---------------------------- REgresión toda la muestra ---------------------------------- #
      
      
      GIRF_data_aux_ps_90 = GIRF_data_aux[,eval(90*19+1):eval(91*19)]
      colnames_salida_state = colnames(GIRF_data_aux_ps_90)
      salida_state_ps_90 = NULL
      for(i in 1:ncol(GIRF_data_aux_ps_90)){
         model_aux = summary(lm(GIRF_data_aux_ps_90[,i]~ 1+state_variable_lm))
         model_aux
         beta = model_aux$coefficients[2,1]
         std = model_aux$coefficients[2,2]
         pval = model_aux$coefficients[2,4]
         salida_aux = cbind(beta,std, pval)
         salida_state_ps_90 = rbind(salida_state_ps_90, salida_aux)
      }
      rownames(salida_state_ps_90) = colnames_salida_state
      salida_state_ps_90[c(2,5,10,11,12,13,15,17,19),]
      
      
      GIRF_data_aux_ns_10 = GIRF_data_aux[,eval(10*19+1):eval(11*19)]
      colnames_salida_state = colnames(GIRF_data_aux_ns_10)
      salida_state_ns_10 = NULL
      for(i in 1:ncol(GIRF_data_aux_ns_10)){
         model_aux = summary(lm(GIRF_data_aux_ns_10[,i]~ 1+state_variable_lm))
         model_aux
         beta = model_aux$coefficients[2,1]
         std = model_aux$coefficients[2,2]
         pval = model_aux$coefficients[2,4]
         salida_aux = cbind(beta, std, pval)
         salida_state_ns_10 = rbind(salida_state_ns_10, salida_aux)
      }
      rownames(salida_state_ns_10) = colnames_salida_state
      salida_state_ns_10[c(2,5,10,11,12,13,15,17,19),]
      plot(GIRF_data_aux_ns_10[,10], type = 'l')      
      
# ----------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------- #      
           
   GIRF_data_aux_ps_90 = GIRF_data_aux[,eval(90*19+1):eval(91*19)]
   colnames_salida_state = colnames(GIRF_data_aux_ps_90)
   salida_state_ps_90_rolling = NULL
   for(j in 0:eval(nrow(GIRF_data_aux_ps_90)-30)){
   salida_state_ps_90 = NULL
   for(i in 1:ncol(GIRF_data_aux_ps_90)){
      model_aux = summary(lm(GIRF_data_aux_ps_90[eval(1+j):eval(30+j),i]~ 1+state_variable_lm[eval(1+j):eval(30+j)]))
      model_aux
      beta = model_aux$coefficients[2,1]
      std = model_aux$coefficients[2,2]
      pval = model_aux$coefficients[2,4]
      salida_aux = cbind(beta,std, pval)
      salida_state_ps_90 = rbind(salida_state_ps_90, salida_aux)
   }
   rownames(salida_state_ps_90) = colnames_salida_state
   aux = salida_state_ps_90
   salida_state_ps_90_rolling = cbind(salida_state_ps_90_rolling,aux)
   }
   salida_state_ps_90_rolling
   write.csv(salida_state_ps_90_rolling, 'betas_ict_ps_90_rolling.csv')
   
   GIRF_data_aux_ps_75 = GIRF_data_aux[,eval(75*19+1):eval(76*19)]
   colnames_salida_state = colnames(GIRF_data_aux_ps_75)
   salida_state_ps_75_rolling = NULL
   for(j in 0:eval(nrow(GIRF_data_aux_ps_75)-30)){
   salida_state_ps_75 = NULL
   for(i in 1:ncol(GIRF_data_aux_ps_75)){
      model_aux = summary(lm(GIRF_data_aux_ps_75[eval(1+j):eval(30+j),i]~ 1+state_variable_lm[eval(1+j):eval(30+j)]))
      model_aux
      beta = model_aux$coefficients[2,1]
      std = model_aux$coefficients[2,2]
      pval = model_aux$coefficients[2,4]
      salida_aux = cbind(beta,std, pval)
      salida_state_ps_75 = rbind(salida_state_ps_75, salida_aux)
   }
   rownames(salida_state_ps_75) = colnames_salida_state
   aux = salida_state_ps_75
   salida_state_ps_75_rolling = cbind(salida_state_ps_75_rolling, aux)
   }
   write.csv(salida_state_ps_75_rolling, 'betas_ict_ps_75_rolling.csv')
   
   
   GIRF_data_aux_ns_25 = GIRF_data_aux[,eval(25*19+1):eval(26*19)]
   colnames_salida_state = colnames(GIRF_data_aux_ns_25)
   salida_state_ns_25_rolling = NULL
   for(j in 0:eval(nrow(GIRF_data_aux_ns_25)-30)){
   salida_state_ns_25 = NULL
   for(i in 1:ncol(GIRF_data_aux_ns_25)){
      model_aux = summary(lm(GIRF_data_aux_ns_25[eval(1+j):eval(30+j),i]~ 1+state_variable_lm[eval(1+j):eval(30+j)]))
      model_aux
      beta = model_aux$coefficients[2,1]
      std = model_aux$coefficients[2,2]
      pval = model_aux$coefficients[2,4]
      salida_aux = cbind(beta, std, pval)
      salida_state_ns_25 = rbind(salida_state_ns_25, salida_aux)
   }
   rownames(salida_state_ns_25) = colnames_salida_state
   aux = salida_state_ns_25
   salida_state_ns_25_rolling = cbind(salida_state_ns_25_rolling, aux)
   }
   write.csv(salida_state_ns_25_rolling, 'betas_ict_ps_25_rolling.csv')
   
   
   GIRF_data_aux_ns_10 = GIRF_data_aux[,eval(10*19+1):eval(11*19)]
   colnames_salida_state = colnames(GIRF_data_aux_ns_10)
   salida_state_ns_10_rolling = NULL
   for(j in 0:eval(nrow(GIRF_data_aux_ns_10)-30)){
   salida_state_ns_10 = NULL
   for(i in 1:ncol(GIRF_data_aux_ns_10)){
      model_aux = summary(lm(GIRF_data_aux_ns_10[eval(1+j):eval(30+j),i]~ 1+state_variable_lm[eval(1+j):eval(30+j)]))
      model_aux
      beta = model_aux$coefficients[2,1]
      std = model_aux$coefficients[2,2]
      pval = model_aux$coefficients[2,4]
      salida_aux = cbind(beta, std, pval)
      salida_state_ns_10 = rbind(salida_state_ns_10, salida_aux)
   }
   rownames(salida_state_ns_10) = colnames_salida_state
   aux = salida_state_ns_10
   salida_state_ns_10_rolling = cbind(salida_state_ns_10_rolling, aux)
   }
   write.csv(salida_state_ns_10_rolling, 'betas_ict_ps_10_rolling.csv')
   
   
   # --------------- Approximation with Panel Model ------------------------- # 
   
   panel_data = NULL 
   index_dates_panel = index_dates_data_state_final[2:48]
   names(index_dates_panel) = "Date"
   shock_index = seq(1,101,1)
   #Panel Data: 
   for(i in 1:nrow(GIRF_state_outputgap)){
      panel_data_aux = NULL
      for(j in seq(1,ncol(GIRF_state_outputgap),19)){
         aux = GIRF_state_outputgap[i,j:eval(j+18)]
         panel_data_aux = rbind(panel_data_aux, aux)
      }
   panel_data_aux = cbind(index_dates_panel[i],shock_index,panel_data_aux, mp_shock_2,positive_shock_sign,state_variable_lm)
   panel_data = rbind(panel_data, panel_data_aux)   
   }
   panel_data = xts(panel_data, as.Date(panel_data[,1]))
   library(plm)
   library(gplots)
   colnames(panel_data)
   panel_data = pdata.frame(panel_data, index = c("shock_index", "V1"))
   table(index(panel_data), useNA = "ifany")
   plotmeans(CSGIRF_int_t_M~V1, data = panel_data)
   
   formula_panel = CSGIRF_inf_t_M ~1 + state_variable_lm + positive_shock_sign
   
   pooltest(formula_panel, data = panel_data, model = "within")
   pFtest(formula_panel, data = panel_data, model = "within")
   plmtest(formula_panel, data = panel_data, effect =  "individual", type = "bp", index = "V1")
   plmtest(formula_panel, data = panel_data, effect = "time", type = "bp")
   plmtest(formula_panel, data = panel_data, effect = "twoways", type = "bp")
   
   phtest(formula_panel,data= panel_data, model = c("random","within"), index = "V1")
   summary(plm(formula_panel,
               data= panel_data, model = "random"))
   summary(plm(formula_panel,
               data= panel_data, model = "pooling" ))
   summary(plm(formula_panel,
               data= panel_data, model = "within", effect = "individual"))
   summary(plm(formula_panel,
              data= panel_data, model = "within", effect = "time"))
   
   summary(plm(CSGIRF_inf_t_M ~state_variable_lm,
               data= panel_data, model = "fd", effect = "individual"))
   
   
   histogram(as.vector(GIRF_int_t[,20]))
   length(state_var_inf_girf)-20
   library(midasr)
   midas_r()