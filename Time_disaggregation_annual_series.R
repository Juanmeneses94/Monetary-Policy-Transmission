# --------------------------------------------- #
# Time Disaggregation of Annual series
# Juan Camilo Meneses Cortés 
# Updated: 27-01-2021
# --------------------------------------------- #

library(readxl)
library(tempdisagg)


#Working Directory
setwd("C:/Users/juanm/OneDrive - Universidad de los Andes/Juan Camilo Meneses/Personal/Tesis PEG/Datos")

#Read Data: 
data_average = read_excel("Datos.xlsx", sheet = "Annual Data", range = "A2:AA24", col_names = T, skip = 0)
names(data_average)

data_average = na.omit(data_average)

data_average_ts = ts(data_average, start = c(2002), frequency = 1)

names_average = colnames(data_average_ts)
names_average = names_average[-1]

salida_average = NULL
for(i in 2:ncol(data_average_ts)){
mod11 <- td(data_average_ts[,i] ~ 1, conversion = "average" ,method = "chow-lin-minrss-ecotrim" , to = "quarterly") 
summary(mod11)  # summary statistics 
#plot(predict(mod11))  
aux  = predict(mod11)
plot(aux)
salida_average = cbind(salida_average,aux)
}
colnames(salida_average) = names_average
write.csv(salida_average, "trimestralización_promedio_var_anuales.csv")

salida_suma = NULL
for(i in 2:ncol(data_average_ts)){
  mod11 <- td(data_average_ts[,i] ~ 1, conversion = "sum" ,method = "chow-lin-minrss-ecotrim" , to = "quarterly") 
  summary(mod11)  # summary statistics 
  #plot(predict(mod11))  
  aux  = predict(mod11)
  plot(aux)
  salida_suma = cbind(salida_suma,aux)
}
colnames(salida_suma) = names_average
write.csv(salida_suma, "trimestralización_suma_var_anuales.csv")


#Informality
informality = read_excel("Datos.xlsx", sheet = "Annual Data", range = "AB6:AB11", col_names = F, skip = 0)
informality = ts(informality, start = c(2001), frequency = 1)
names(informality) = c("informality")
mod_inf <- td(informality ~ 1, conversion = "average" ,method = "chow-lin-minrss-ecotrim" , to = "quarterly") 
aux  = predict(mod_inf)
plot(aux)
write.csv(aux, "trimestralización_informalidad.csv")






