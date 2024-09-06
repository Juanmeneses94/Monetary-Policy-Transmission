# Monetary-Policy-Transmission
R codes and csv/xlsx databases that support econometric/statistical methdologies used in the research paper: "Index of Conditions of Monetary Policy Transmission in Colombia"
Ms Economics. Univerisdad de los Andes.
June 2021
Author: Juan Camilo Meneses

The methodologies were used to estimate non-linear impulse-response functions over a three-equational macroeconomic model of the Colombian economy, with:

1. A demand function "IS", that describes output gap dynamics.
2. A Phillips-Curve, that describes inflation dynamics.
3. A monetary policy "rule" a la Taylor.    

Data includes +164 macroeconomic/financial variables of the Colombian economy including national statistics, survey data, market data, and transformation of those variables (such as annual growht rates, moving averages, conversion to dollar amounts, and so on). A dictionary is found in sheet "Variables" in "Datos.xlsx".

The following are the main codes and steps used in the estimation:

# Output_gap_estimates.R: 
Code main output is the estimation of the output gap of the Colombian economy, using various methdologies common in the economic literature. It includes:

Pre-processing: 
1. Seasonal Adjustment of time series using "seasonal" library. 
2. H-P filters using "mFilter" library.

Output Gap: 
Using GDP series, univariate filters:
1. H-P Filter with "mFilter" library, including filter over de trend-cycle component.
2. Extended H-P Filter, which includes function to "extend" time series using forecasts and backcasts.
3. Never Use H-P Filter, by Hamilton (2017), using library "neverhpfilter".
Using structural approaches:
1. Cobb-Douglas production functions: includes estimates of NAIRU using splines, and trend-cycle decomposition of unemployment series. Also includes Kalman Filter for capaicty utilization series to obatin NAICU.
2. Kalman filter using "dse" and "DLM" libraries. State-space equation include output gap, core inflation, real exchange rate depreciation, real monetary policy rate.

# Time_disaggregation_time_series.R:
Part of the pre-processing of data before estimates. It uses library "tempdisagg" on annual data to transform to quarterly data for several variables, including World Bank indexes, Gini index, world gdp, among others. 

# Nonlinearity, state-dependency and asymmetry in the MTM - GIRF
This is the core code for the estimation of the three-equational model of the Colombian economy.
As a first step, it includes a block for unit root test under three methodologies (Dickey-Fuller, ERS, and KPSS, together with Enders & Ludlow (2000).

Methodologies included in the code are: 
1. Non Linear Least Squares Estimates using libraries "nlme" and "nlstools".
2. Testing for Nonlinearity and State Dependency (LM Type-Tests). 
3. Selection of transition functions (STAR, LSTAR1 or LSTAR 2).
4. Estimates of STR models and diagnosis tests using parallel processing. 
5. Estimates of GIRF (densities). The last step is computed following Koop (1996), with the following algorithm:

Koop (1996) for multivariate models:

1. Get the variance covariance matrix of a MMT model residuals (K*T) ~ (3*68)
2. Compute Cholesky factorization of empirical resid var-cov, its inverse and orthogonalize residuals using Cholesky factor
3. Draw randomly from the independent (K*T) resid matrix the interest rate shock (K - dimensional), and use Cholesky factor to return to dependence. Poztek uses percentile 0 to 100 in this step of the empirical distribution of interest rate equation. Here the propose approach is to use the decomposition, take the percentiles of the interes rate equation and then return to dependence. 
4. Select shock (percentile of interest equation distribution derived in 4) and a point in time for the model variable (from t+3.....T).
5. For a given horizon N (here set at 16 periods -4 years-), randomly sample with (N+1)*R (R = repetitions/10000) from the K*T innovation matrix 
6. Use first N random innovations to iterate on nonlinear equations to get y(vt, wt-1), and N+1 random innovations to iterate on nonlinear model to obtain y(wt-1) 
7. Take the difference of y(vt, wt-1) - y(wt-1) for each N. (this is the GIRF)
8. Standardize GIRF by initial shock (chosen in step 4). 
9. Repeat step 6-8 R times to form the distribution of Monte Carlo simulation. 
10.Compute measures of interest (i.e. median, 5th, 20th, 75th, 80th pctiles)

