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

# Time_desaggregation_time_series:
Part of the pre-processing of data before estimates. It uses library "tempdisagg" on annual data to transform to quarterly data for several variables, including World Bank indexes, Gini index, world gdp, among others. 

# 

