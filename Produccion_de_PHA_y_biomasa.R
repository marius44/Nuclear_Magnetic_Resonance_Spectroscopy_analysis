# Plots crecimento y PHAs Andrea
# Junio 2025
# Andy datos UAMI

library("ggplot2")
library("dplyr")
library("tidyr")
library("patchwork") # Para ver 2 ejes Y
library("ggalt")

setwd(dir = "D:/Proyectos_y_revisiones/Académicos/Andrea_maestria/plots_crecimiento/")
getwd()

# Importar datos
cinetica <- read.csv(file = "Datos_para_plots.csv", header = T)
str(cinetica)
View(cinetica)

# Remobra las columnas
names(cinetica) <- gsub("\\.+", "_", names(cinetica))  # reemplaza múltiples puntos por "_"
names(cinetica) <- gsub("_$", "", names(cinetica))     # elimina guiones bajos al final

#View(cinetica)

# UNa solo grafica

crecimiento_80g <- ggplot(cinetica, 
                          aes(x=Tiempo_hr, y=OD_600_80_g_L_NaCl))+
  geom_point() +
  geom_line()+
  ylab(label = "OD 80/L NaCl")+
  xlab(label = "Tiempo (h)")
crecimiento_80g

###LONGER
longer <- cinetica %>%
  pivot_longer(cols = -'Tiempo_hr',
               names_to = "variable",
               values_to = "valor"
               )
View(longer)
str(longer)
#FUNCION PARA SUAVIZAR Y REDONDEAR LINEAS
# https://stackoverflow.com/questions/43771900/use-curved-lines-in-bumps-chart
ggpchip = function(formula, data, weights) structure(pracma::pchipfun(data$x, data$y), class='ggpchip')
predict.ggpchip = function(object, newdata, se.fit=F, ...) {
  fit = unclass(object)(newdata$x)
  if (se.fit) list(fit=data.frame(fit, lwr=fit, upr=fit), se.fit=fit * 0) else fit
}

# SOLO DO
long_od <- longer %>% 
  filter(grepl("^OD", variable)) %>%
  ggplot(aes(x = Tiempo_hr, y = valor, color = variable)) +
  geom_point() +
  #geom_line() +
  geom_smooth(method='ggpchip', se=F, size = 1)+
  labs(title = "Crecimiento", x = "Tiempo (h)", y = "Valor medido", color = "") +
  ylab(label = "OD (600 nm)")+
  xlab(label = "Tiempo (h)")+
  theme_minimal()
long_od
?geom_smooth
# SOLO PHA
long_pha <- longer %>%
  filter(grepl("^PHA_", variable)) %>%     # filtra variables que empiecen por 'pha_'
  ggplot(aes(Tiempo_hr, valor,
             color = variable, group = variable)) +
  geom_point() +
  #geom_line() +
  geom_smooth(method='ggpchip', se=F, size = 1)+
  labs(title = "Producción PHAs", x = "Tiempo (h)", y = "Valor medido", color = "") +
  ylab(label = "PHA (mg)")+
  xlab(label = "Tiempo (h)")+
  theme_minimal()
long_pha

# SOLO PH
long_ph <- longer %>%
  filter(grepl("^pH_", variable)) %>%     # filtra variables que empiecen por 'pha_'
  ggplot(aes(Tiempo_hr, valor,
             color = variable, group = variable)) +
  geom_point() +
  #geom_line() +
  geom_smooth(method='ggpchip', se=F, size = 1)+
  labs(title = "pH", x = "Tiempo (h)", y = "Valor medido", color = "") +
  ylab(label = "pH")+
  xlab(label = "Tiempo (h)")+
  theme_minimal()
long_ph

#MASA
long_msc <- longer %>%
  filter(grepl("^MSC_", variable)) %>%     # filtra variables que empiecen por 'pha_'
  ggplot(aes(Tiempo_hr, valor,
             color = variable, group = variable)) +
  geom_point() +
  #geom_line() +
  geom_smooth(method='ggpchip', se=F)+
  labs(title = "Masa seca celular", x = "Tiempo (h)", y = "Valor medido", color = "") +
  ylab(label = "MSC (mg)")+
  xlab(label = "Tiempo (h)")+
  theme_minimal()
long_msc

################################################################################
# Filtrar datos 80 g y %PHA
msc <- longer %>% filter(variable == "MSC_80_g_L_NaCl")
pha <- longer %>% filter(variable == "REND_PHA_MSC_80_g_L_NaCl")


# Verificar
if (nrow(msc) == 0 || nrow(pha) == 0) {
  stop("MSC_80_g_L_NaCl o PHA_80_g_L_NaCl no tienen datos.")
}

# Escalar
max_msc <- max(msc$valor, na.rm = TRUE)
max_pha <- max(pha$valor, na.rm = TRUE)
escala <- max_msc / max_pha

# Transformar y etiquetar
msc <- msc %>% mutate(valor_esc = valor, variable_plot = "MSC")
pha <- pha %>% mutate(valor_esc = valor * escala, variable_plot = "PHA")

# Unir y promediar por tiempo para evitar duplicados
df_plot <- bind_rows(msc, pha) %>%
  group_by(variable_plot, Tiempo_hr) %>%
  summarize(valor_esc = mean(valor_esc, na.rm = TRUE), .groups = "drop") %>%
  arrange(variable_plot, Tiempo_hr)

# Graficar
ggplot(df_plot, aes(x = Tiempo_hr, y = valor_esc, color = variable_plot)) +
  geom_point() +
  geom_smooth(method = "ggpchip", se = FALSE) +
  scale_y_continuous(
    name = "MSC (mg)",
    sec.axis = sec_axis(~ . / escala, name = "PHA (% MSC)")
  ) +
  labs(x = "Tiempo (h)", title = "Biomasa contra producción (%) de PHA  (80 g/L NaCl)", color = "") +
  theme_minimal()



# Filtrar datos 100 g y %PHA
msc <- longer %>% filter(variable == "MSC_100_g_L_NaCl")
pha <- longer %>% filter(variable == "REND_PHA_MSC_100_g_L_NaCl")


# Verificar
if (nrow(msc) == 0 || nrow(pha) == 0) {
  stop("MSC_100_g_L_NaCl o REND_PHA_MSC_100_g_L_NaCl no tienen datos.")
}

# Escalar
max_msc <- max(msc$valor, na.rm = TRUE)
max_pha <- max(pha$valor, na.rm = TRUE)
escala <- max_msc / max_pha

# Transformar y etiquetar
msc <- msc %>% mutate(valor_esc = valor, variable_plot = "MSC")
pha <- pha %>% mutate(valor_esc = valor * escala, variable_plot = "PHA")

# Unir y promediar por tiempo para evitar duplicados
df_plot <- bind_rows(msc, pha) %>%
  group_by(variable_plot, Tiempo_hr) %>%
  summarize(valor_esc = mean(valor_esc, na.rm = TRUE), .groups = "drop") %>%
  arrange(variable_plot, Tiempo_hr)

# Graficar
ggplot(df_plot, aes(x = Tiempo_hr, y = valor_esc, color = variable_plot)) +
  geom_point() +
  geom_smooth(method = "ggpchip", se = FALSE) +
  scale_y_continuous(
    name = "MSC (mg)",
    sec.axis = sec_axis(~ . / escala, name = "PHA (% MSC)")
  ) +
  labs(x = "Tiempo (h)", title = "Biomasa contra producción (%) de PHA  (100 g/L NaCl)", color = "") +
  theme_minimal()


# Filtrar datos 120 g y %PHA
msc <- longer %>% filter(variable == "MSC_120_g_L_NaCl")
pha <- longer %>% filter(variable == "REND_PHA_MSC_120_g_L_NaCl")


# Verificar
if (nrow(msc) == 0 || nrow(pha) == 0) {
  stop("MSC_120_g_L_NaCl o REND_PHA_120_g_L_NaCl no tienen datos.")
}

# Escalar
max_msc <- max(msc$valor, na.rm = TRUE)
max_pha <- max(pha$valor, na.rm = TRUE)
escala <- max_msc / max_pha

# Transformar y etiquetar
msc <- msc %>% mutate(valor_esc = valor, variable_plot = "MSC")
pha <- pha %>% mutate(valor_esc = valor * escala, variable_plot = "PHA")

# Unir y promediar por tiempo para evitar duplicados
df_plot <- bind_rows(msc, pha) %>%
  group_by(variable_plot, Tiempo_hr) %>%
  summarize(valor_esc = mean(valor_esc, na.rm = TRUE), .groups = "drop") %>%
  arrange(variable_plot, Tiempo_hr)

# Graficar
ggplot(df_plot, aes(x = Tiempo_hr, y = valor_esc, color = variable_plot)) +
  geom_point() +
  geom_smooth(method = "ggpchip", se = FALSE) +
  scale_y_continuous(
    name = "MSC (mg)",
    sec.axis = sec_axis(~ . / escala, name = "PHA (% MSC)")
  ) +
  labs(x = "Tiempo (h)", title = "Biomasa contra producción (%) de PHA  (120 g/L NaCl)", color = "") +
  theme_minimal()

################################################################################
# AHORA TODO
# Filtrar MSC y PHA
msc <- longer %>% filter(variable %in% c("MSC_80_g_L_NaCl", "MSC_100_g_L_NaCl", "MSC_120_g_L_NaCl"))
rend_pha <- longer %>% filter(variable %in% c("REND_PHA_MSC_80_g_L_NaCl", "REND_PHA_MSC_100_g_L_NaCl", "REND_PHA_MSC_120_g_L_NaCl"))

# Calcular escala con máximos reales
max_msc <- max(msc$valor, na.rm = TRUE)
rend_max_pha <- max(rend_pha$valor, na.rm = TRUE)
escala <- max_msc / rend_max_pha

# Escalar PHA para que comparta eje con MSC
rend_pha <- rend_pha %>%
  mutate(valor_esc = valor * escala)

# MSC ya está en la escala deseada
msc <- msc %>%
  mutate(valor_esc = valor)

# Unir los datos
df_plot <- bind_rows(msc, rend_pha) %>%
  mutate(grupo = case_when(
    grepl("MSC", variable) ~ gsub("_g_L_NaCl", "", variable),
    grepl("REND_PHA", variable) ~ gsub("_g_L_NaCl", "", variable)
  ),
  tipo = ifelse(grepl("MSC", variable), "MSC", "REND_PHA"))

# Graficar
ggplot(df_plot, aes(x = Tiempo_hr, y = valor_esc, color = grupo)) +
  geom_point() +
  geom_smooth(method = "ggpchip", se = FALSE) +
  scale_y_continuous(
    name = "MSC (mg)",
    sec.axis = sec_axis(~ . / escala, name = "PHA (% MSC)")
  ) +
  labs(x = "Tiempo (h)", title = "Producción y biomasa", color = " ") +
  theme_minimal()


