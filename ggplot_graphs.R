
# Cargamos los datos
library(ggplot2)

################
# Gráfico base #
################

load(file = "josep_Arriaga.Rdata")

#datos_arriaga <- read.csv(file = "datos/datos_arriaga.csv", sep = ";", header = T)

datos_arriaga <- rbind.data.frame(cbind.data.frame(age = 0:100,
                                                   delta = cc_2010_mean,
                                                   year = 2010),
                                  cbind.data.frame(age = 0:100,
                                                   delta = cc_2015_mean,
                                                   year = 2015),
                                  cbind.data.frame(age = 0:100,
                                                   delta = cc_2020_mean,
                                                   year = 2020))

png("Figuras/Figura_arriaga.png", height = 5.44, width = 8,
    units = 'in', res = 1500)

ggplot(datos_arriaga, aes(x = age, y = delta, linetype = as.factor(year))) +
  geom_line(size = 0.8) +
  theme_minimal() +
  theme(legend.position = c(0.2, 0.8),
        plot.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(linetype = "Year") +
  xlab("Age") +
  ylab(expression(paste(""[n], Delta[x], sep ="")))

dev.off()

#########################
# Gráfico sensibilities #
#########################

datos_sensi <- read.csv(file = "datos/datos_sensibilidad.csv", sep = ";", header = T)

tiff("Figuras/Figura_sensibilities.tiff", height = 5.44, width = 8,
     units = 'in', compression = "lzw", res = 1500)

ggplot(datos_sensi, aes(x = age, y = s)) +
  geom_point(shape = 1) +
  geom_line(aes(x = age, y = s1), col = "#800000")+
  geom_line(aes(x = age, y = s2), col = "#0000FF") +
  geom_line(aes(x = age, y = s1_s2), linetype = "dashed") +
  geom_line(aes(x = age, y = s1_d), col = "#800000", alpha = 0.3, linetype = "dashed") +
  geom_line(aes(x = age, y = s2_d), col = "#0000FF", alpha = 0.3, linetype = "dashed") +
  xlab("Age") +
  ylab(expression(paste(""[n], Delta[x], " and sensibilities", sep =""))) +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

dev.off()

##############################
# Gráfico Pirámide Población #
##############################
library(dplyr)

load(file = "US_Exp.Rdata")

datos_piramide_bis <-  sp_E %>%
  filter(Year == 2003)

# Sumar las variables para obtener valores positivos y negativos
datos_piramide_bis$Male <- abs(datos_piramide_bis$Male)
datos_piramide_bis$Female <- -abs(datos_piramide_bis$Female)

tiff("Figuras/Poblacion USA.tiff", height = 5.44, width = 8,
     units = 'in', compression = "lzw", res = 1500)

ggplot(datos_piramide_bis, aes(x=Age)) +
  geom_bar(aes(y=Male), fill="#0000FF", stat="identity", width=0.75, position="identity") +
  geom_bar(aes(y=Female), fill = "#800000", stat="identity", width=0.75, position="identity") +
  coord_flip() +
  scale_y_continuous(breaks = seq(-2400000, 2400000, 800000),
                     labels = paste0(as.character(c(seq(2400, 0, -800), seq(400, 2400, 800))), "m"),
                     limits = c(-2400000, 2400000)) +
  xlim(c(0,100)) +
  xlab("Age") + ylab("Population number of U.S.A") +
  theme(plot.title = element_text(hjust = 0.5,size=12,face="bold"),
        panel.border = element_rect(colour = "black",fill=NA,size=1),
        axis.text = element_text(size=12,face="bold"),
        panel.background = element_rect(fill= NA)) +
  annotate("text", x = 80, y = 1400000, label = "Male", hjust = 0, col = "#0000FF") +
  annotate("text", x = 80, y = -1800000, label = "Female", hjust = 0, col = "#800000")

dev.off()

#######################
# Gráfico Defunciones #
#######################

defunciones <- read.csv(file = "datos/Defunciones.csv", sep = ";", header = T)

defunciones <- defunciones %>% filter(year_id %in% c(2010, 2015, 2020))

defunciones$gender <- factor(defunciones$gender, levels = c("Male", "Female"))


tiff("Figuras/Defunciones.tiff", height = 5.44, width = 8,
     units = 'in', compression = "lzw", res = 1500)

ggplot(defunciones, aes(x = age_id, y = Dxt_HMD, color = as.factor(year_id))) +
  geom_line() +
  stat_smooth(geom = "line", method = "loess",
              alpha = 0.2, size = 0.5, se = F, span = 0.50) +
  facet_wrap("gender") +
  scale_color_manual(values=c("orange4", "#293352", "olivedrab4")) +
  ylab("Number of deaths") + xlab("Age") +
  labs(color = "Year") +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

dev.off()


#################################
#  Gráfico mx y descomposición  #
#################################

load(file = "Josep.Rdata")


library(tidyr)
df_M <- as.data.frame(pivot_longer(US_male,
                                   cols = mx:"mx-V01-Y89",
                                   names_to = "cause",
                                   values_to = "mx_cause"))
df_F <- pivot_longer(US_fema,
                     cols = mx:"mx-V01-Y89",
                     names_to = "cause",
                     values_to = "mx_cause")

df_M_1 <- df_M[df_M$cause %in% c("mx", "mx-A00-B99", "mx-C00-D48",
                                 "mx-D50-D89", "mx-E00-E88") &
                 df_M$Per == 2015, ]

df_M_1 <- df_M_1[df_M_1$mx_cause > 0, ]

ggplot(df_M_1, aes(x = Age, y = log(mx_cause), color = as.factor(cause))) +
  geom_point()




#################################
#  descomposición de Arriaga    #
#  de las diferentes causas     #
#  por grupos de edad           #
#################################

#load(file = "Grafico_Causes_Decomp.RData")

#write.table(arriaga_causes_decomp, file = "ayuda.csv", sep = ";", row.names = F)

df <- read.csv(file = "ayuda.csv", sep = ";", header = T)

df$age <- 0:100

library(tidyr)
library(ggplot2)
library(ggridges)
df <- pivot_longer(df,
                   cols = A00.B99:V01.Y89,
                   names_to = "tipo",
                   values_to = "total")

df <- as.data.frame(df)

df$Total_n <- floor(df$total * 100000)

grupo_acum <- edad_acum <- c()
for (k in 1:dim(df)[1]){
  
  edad_rep <- rep(df$edad[k], abs(df$Total_n[k]))
  if (df$Total_n[k] < 0) {
    edad_rep <- edad_rep * (-1)
  }
  causa_rep <- rep(df$tipo[k], abs(df$Total_n[k]))
  edad_acum <- c(edad_acum, edad_rep)
  grupo_acum <- c(grupo_acum, causa_rep)
  
}


ggplot(df, aes(x = age, y = - total)) +
  geom_point() +
  facet_wrap("tipo") +
  geom_smooth()







