library(VGAM)
library(VGAMdata)
library(MVN)


data("xs.nz")

variables_modelo <- c(
  "sbp",
  "dbp",
  "smokenow",
  "sex",
  "weight",
  "height",
  "age"
)

base_xs <- xs.nz[
  complete.cases(xs.nz[, variables_modelo]),
  variables_modelo
]

base_xs$tratamiento <- factor(base_xs$smokenow)
base_xs$sex <- factor(base_xs$sex)

familia_binormal <- VGAM::binormal(
  lmean1 = "identitylink",
  lmean2 = "identitylink",
  lsd1   = "loglink",
  lsd2   = "loglink",
  lrho   = "rhobitlink"
)

modelo_parametrico <- VGAM::vgam(
  cbind(sbp, dbp) ~
    tratamiento +
    sex +
    weight +
    height +
    age,
  family = familia_binormal,
  data = base_xs,
  trace = TRUE
)

modelo_semiparametrico <- VGAM::vgam(
  cbind(sbp, dbp) ~
    tratamiento +
    sex +
    weight +
    height +
    VGAM::sm.ps(age),
  family = familia_binormal,
  data = base_xs,
  trace = TRUE
)



comparacion_modelos <- data.frame(
  modelo = c(
    "Paramétrico",
    "Semiparamétrico"
  ),
  logLik = c(
    as.numeric(logLik(modelo_parametrico)),
    as.numeric(logLik(modelo_semiparametrico))
  ),
  AIC = c(
    AIC(modelo_parametrico),
    AIC(modelo_semiparametrico)
  ),
  BIC = c(
    BIC(modelo_parametrico),
    BIC(modelo_semiparametrico)
  ),
  n = rep(nrow(base_xs), 2)
)

comparacion_modelos$delta_AIC <-
  comparacion_modelos$AIC -
  min(comparacion_modelos$AIC)

comparacion_modelos$delta_BIC <-
  comparacion_modelos$BIC -
  min(comparacion_modelos$BIC)

print(
  comparacion_modelos,
  digits = 10,
  row.names = FALSE
)

plot(
  modelo_semiparametrico,
  se = TRUE,
  residuals = FALSE,
  rugplot = TRUE
)














edad <- seq(
  min(base_xs$age),
  max(base_xs$age),
  length.out = 200
)

nuevo <- data.frame(
  tratamiento = factor(
    levels(base_xs$tratamiento)[1],
    levels = levels(base_xs$tratamiento)
  ),
  sex = factor(
    levels(base_xs$sex)[1],
    levels = levels(base_xs$sex)
  ),
  weight = mean(base_xs$weight),
  height = mean(base_xs$height),
  age = edad
)

terminos <- predict(
  modelo_semiparametrico,
  newdata = nuevo,
  type = "terms"
)

colnames(terminos)

idx <- grep(
  "sm.ps",
  colnames(terminos),
  fixed = TRUE
)

funcion_suave <- data.frame(
  age = edad,
  f_sbp = terminos[, idx[1]],
  f_dbp = terminos[, idx[2]]
)

f_sbp <- splinefun(
  funcion_suave$age,
  funcion_suave$f_sbp
)

f_dbp <- splinefun(
  funcion_suave$age,
  funcion_suave$f_dbp
)

funcion_suave


# Construir funciones cúbicas por tramos
f_sbp <- splinefun(
  funcion_suave$age,
  funcion_suave$f_sbp,
  method = "fmm"
)

f_dbp <- splinefun(
  funcion_suave$age,
  funcion_suave$f_dbp,
  method = "fmm"
)

# Extraer coeficientes polinómicos de cada tramo
extraer_ecuaciones <- function(f) {
  
  z <- environment(f)$z
  
  data.frame(
    inferior = z$x[-length(z$x)],
    superior = z$x[-1],
    b0 = z$y[-length(z$y)],
    b1 = z$b[-length(z$b)],
    b2 = z$c[-length(z$c)],
    b3 = z$d[-length(z$d)]
  )
}

ecuacion_sbp <- extraer_ecuaciones(f_sbp)
ecuacion_dbp <- extraer_ecuaciones(f_dbp)

imprimir_ecuaciones <- function(tabla, nombre) {
  
  for (i in seq_len(nrow(tabla))) {
    
    cat(
      sprintf(
        "%s(age) = %.8f %+.8f*(age - %.5f) %+.8f*(age - %.5f)^2 %+.8f*(age - %.5f)^3,\n",
        nombre,
        tabla$b0[i],
        tabla$b1[i],
        tabla$inferior[i],
        tabla$b2[i],
        tabla$inferior[i],
        tabla$b3[i],
        tabla$inferior[i]
      )
    )
    
    cat(
      sprintf(
        "para %.5f <= age < %.5f\n\n",
        tabla$inferior[i],
        tabla$superior[i]
      )
    )
  }
}

imprimir_ecuaciones(ecuacion_sbp, "f_SBP")
imprimir_ecuaciones(ecuacion_dbp, "f_DBP")



# ---------------------------------------------------------
# Normalidad bivariada de SBP y DBP
# ---------------------------------------------------------

library(MVN)

datos_normalidad <- base_xs[, c("sbp", "dbp")]

# ---------------------------------------------------------
# Henze-Zirkler
# ---------------------------------------------------------

normalidad_hz <- mvn(
  datos_normalidad,
  mvn_test = "hz"
)

normalidad_hz


# ---------------------------------------------------------
# Mardia
# ---------------------------------------------------------

normalidad_mardia <- mvn(
  datos_normalidad,
  mvn_test = "mardia"
)

normalidad_mardia
