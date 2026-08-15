# ============================================================
# Paquetes
# ============================================================

library(VGAM)
library(VGAMdata)
library(vgamCausal)
library(MatchIt)
library(cobalt)
library(ggplot2)
library(boot)


# ============================================================
# Semilla
# ============================================================

set.seed(1987)


# ============================================================
# Datos
# ============================================================

data("xs.nz")


# ============================================================
# Selección de variables
# ============================================================

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


# ============================================================
# Tratamiento
# ============================================================

base_xs$tratamiento <- factor(
  base_xs$smokenow,
  levels = c(0, 1),
  labels = c(
    "No_fumador",
    "Fumador"
  )
)


base_xs$sex <- factor(base_xs$sex)


# ============================================================
# Familia VGAM
# ============================================================

familia_binormal <- VGAM::binormal(
  lmean1 = "identitylink",
  lmean2 = "identitylink",
  lsd1   = "loglink",
  lsd2   = "loglink",
  lrho   = "rhobitlink"
)


# ============================================================
# PSM dirigido al ATE
# ============================================================

psm_ate <- vgamc_match_psm(
  data = base_xs,
  treatment = "tratamiento",
  covariates = c(
    "sex",
    "weight",
    "height",
    "age"
  ),
  reference = "No_fumador",
  levels = "Fumador",
  method = "subclass",
  distance = "glm",
  subclass = 5,
  ratio = NULL,
  replace = NULL,
  estimand = "ATE"
)


# ============================================================
# Datos después de la subclasificación
# ============================================================

base_emparejada <- psm_ate$matches[["Fumador"]]$data


# Distribución por subclase y tratamiento

table(
  base_emparejada$subclass,
  base_emparejada$tratamiento
)


# ============================================================
# Balance de covariables
# ============================================================

grafico_balance <- love.plot(
  tratamiento ~
    sex +
    weight +
    height +
    age,
  data = base_emparejada,
  weights = base_emparejada$weights,
  stats = "mean.diffs",
  estimand = "ATE",
  s.d.denom = "pooled",
  binary = "std",
  abs = TRUE,
  thresholds = c(m = 0.10),
  sample.names = c(
    "Before",
    "After"
  ),
  var.names = c(
    sex = "Sexo",
    weight = "Peso",
    height = "Talla",
    age = "Edad"
  ),
  line = TRUE,
  stars = "raw",
  grid = TRUE,
  title = ""
)


# ============================================================
# Guardar gráfico de balance
# ============================================================

ruta_salida <- "C:/RIBG/Causalidad/xs.nz"


if (!dir.exists(ruta_salida)) {
  dir.create(
    ruta_salida,
    recursive = TRUE
  )
}


ggsave(
  filename = file.path(
    ruta_salida,
    "love_plot_ATE.pdf"
  ),
  plot = grafico_balance,
  device = cairo_pdf,
  width = 8,
  height = 5.5,
  units = "in"
)


# ============================================================
# Modelo VGAM
# ============================================================

modelo_psm_vgam_spline <- vgamc_fit_psm_vgam(
  cbind(sbp, dbp) ~
    tratamiento +
    sex +
    weight +
    height +
    VGAM::sm.ps(age),
  matches = psm_ate,
  family = familia_binormal,
  use_weights = TRUE,
  trace = FALSE
)


# ============================================================
# ATE mediante subclases
# ============================================================

ate_spline_subclases <- vgamc_psm_effects(
  fit = modelo_psm_vgam_spline,
  estimands = "ate",
  group = "subclass",
  use_weights = TRUE,
  strict_estimand = TRUE
)


# ============================================================
# Estimaciones puntuales
# ============================================================

resultados_ate <- ate_spline_subclases$effects


resultados_ate_impresion <- resultados_ate


resultados_ate_impresion$estimate <- round(
  resultados_ate_impresion$estimate,
  digits = 3
)


print(
  resultados_ate_impresion,
  row.names = FALSE
)


# ============================================================
# Función para el bootstrap
#
# En cada réplica:
#
# 1. Se remuestrean las observaciones.
# 2. Se vuelve a estimar el propensity score.
# 3. Se vuelven a construir las 5 subclases.
# 4. Se vuelven a calcular los pesos.
# 5. Se vuelve a ajustar el VGAM.
# 6. Se vuelve a estimar el ATE.
# ============================================================

estima_ate_boot <- function(datos, indices) {
  
  # ----------------------------------------------------------
  # Muestra bootstrap
  # ----------------------------------------------------------
  
  datos_b <- datos[
    indices,
    ,
    drop = FALSE
  ]
  
  
  # Evitar nombres de filas duplicados
  
  rownames(datos_b) <- NULL
  
  
  # Mantener niveles originales de los factores
  
  datos_b$tratamiento <- factor(
    datos_b$tratamiento,
    levels = c(
      "No_fumador",
      "Fumador"
    )
  )
  
  
  datos_b$sex <- factor(
    datos_b$sex,
    levels = levels(base_xs$sex)
  )
  
  
  # ----------------------------------------------------------
  # Volver a estimar propensity score y subclases
  # ----------------------------------------------------------
  
  psm_b <- vgamc_match_psm(
    data = datos_b,
    treatment = "tratamiento",
    covariates = c(
      "sex",
      "weight",
      "height",
      "age"
    ),
    reference = "No_fumador",
    levels = "Fumador",
    method = "subclass",
    distance = "glm",
    subclass = 5,
    ratio = NULL,
    replace = NULL,
    estimand = "ATE"
  )
  
  
  # ----------------------------------------------------------
  # Volver a ajustar el VGAM
  # ----------------------------------------------------------
  
  modelo_b <- vgamc_fit_psm_vgam(
    cbind(sbp, dbp) ~
      tratamiento +
      sex +
      weight +
      height +
      VGAM::sm.ps(age),
    matches = psm_b,
    family = familia_binormal,
    use_weights = TRUE,
    trace = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Volver a calcular el ATE
  # ----------------------------------------------------------
  
  ate_b <- vgamc_psm_effects(
    fit = modelo_b,
    estimands = "ate",
    group = "subclass",
    use_weights = TRUE,
    strict_estimand = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Extraer resultados
  # ----------------------------------------------------------
  
  efectos_b <- ate_b$effects
  
  
  ate_sbp <- efectos_b$estimate[
    efectos_b$outcome == "sbp"
  ]
  
  
  ate_dbp <- efectos_b$estimate[
    efectos_b$outcome == "dbp"
  ]
  
  
  # ----------------------------------------------------------
  # Devolver estimaciones
  # ----------------------------------------------------------
  
  c(
    sbp = ate_sbp,
    dbp = ate_dbp
  )
}


# ============================================================
# Bootstrap estratificado por tratamiento
# ============================================================

set.seed(1987)


# Número de réplicas bootstrap

B <- 1000L


boot_ate <- boot(
  data = base_xs,
  statistic = estima_ate_boot,
  R = B,
  strata = base_xs$tratamiento
)


# ============================================================
# Resultados básicos del bootstrap
# ============================================================

print(boot_ate)


# ============================================================
# Verificar réplicas
# ============================================================

replicas_validas <- apply(
  boot_ate$t,
  MARGIN = 1,
  FUN = function(x) {
    all(is.finite(x))
  }
)


cat(
  "\nNúmero total de réplicas:",
  B,
  "\n"
)


cat(
  "Réplicas válidas:",
  sum(replicas_validas),
  "\n"
)


cat(
  "Réplicas con problemas:",
  sum(!replicas_validas),
  "\n\n"
)


# boot.ci requiere estimaciones válidas para calcular BCa

if (any(!replicas_validas)) {
  
  stop(
    paste0(
      "Se encontraron ",
      sum(!replicas_validas),
      " réplicas bootstrap sin estimaciones válidas. ",
      "Revise las réplicas antes de calcular el intervalo BCa."
    )
  )
}


# ============================================================
# Distribución bootstrap
# ============================================================

bootstrap_sbp <- boot_ate$t[, 1]

bootstrap_dbp <- boot_ate$t[, 2]


# ============================================================
# Error estándar bootstrap
# ============================================================

se_boot_sbp <- sd(
  bootstrap_sbp
)


se_boot_dbp <- sd(
  bootstrap_dbp
)


# ============================================================
# Sesgo bootstrap
# ============================================================

bias_boot_sbp <- mean(
  bootstrap_sbp
) - boot_ate$t0[1]


bias_boot_dbp <- mean(
  bootstrap_dbp
) - boot_ate$t0[2]


# ============================================================
# Intervalos de confianza para SBP
#
# Se calculan:
#
# - Percentil
# - BCa
# ============================================================

ic_sbp <- boot.ci(
  boot_ate,
  conf = 0.95,
  type = c(
    "perc",
    "bca"
  ),
  index = 1
)


# ============================================================
# Intervalos de confianza para DBP
# ============================================================

ic_dbp <- boot.ci(
  boot_ate,
  conf = 0.95,
  type = c(
    "perc",
    "bca"
  ),
  index = 2
)


# ============================================================
# Mostrar intervalos completos
# ============================================================

print(ic_sbp)

print(ic_dbp)


# ============================================================
# Extraer intervalo percentil
# ============================================================

ic_percentil_sbp <- c(
  inferior = ic_sbp$percent[4],
  superior = ic_sbp$percent[5]
)


ic_percentil_dbp <- c(
  inferior = ic_dbp$percent[4],
  superior = ic_dbp$percent[5]
)


# ============================================================
# Extraer intervalo BCa
# ============================================================

ic_bca_sbp <- c(
  inferior = ic_sbp$bca[4],
  superior = ic_sbp$bca[5]
)


ic_bca_dbp <- c(
  inferior = ic_dbp$bca[4],
  superior = ic_dbp$bca[5]
)


# ============================================================
# Tabla final
# ============================================================

resultados_ic <- data.frame(
  
  outcome = c(
    "sbp",
    "dbp"
  ),
  
  ATE = c(
    boot_ate$t0[1],
    boot_ate$t0[2]
  ),
  
  SE_bootstrap = c(
    se_boot_sbp,
    se_boot_dbp
  ),
  
  Bias_bootstrap = c(
    bias_boot_sbp,
    bias_boot_dbp
  ),
  
  Percentil_95_inferior = c(
    ic_percentil_sbp["inferior"],
    ic_percentil_dbp["inferior"]
  ),
  
  Percentil_95_superior = c(
    ic_percentil_sbp["superior"],
    ic_percentil_dbp["superior"]
  ),
  
  BCa_95_inferior = c(
    ic_bca_sbp["inferior"],
    ic_bca_dbp["inferior"]
  ),
  
  BCa_95_superior = c(
    ic_bca_sbp["superior"],
    ic_bca_dbp["superior"]
  ),
  
  row.names = NULL
)


# ============================================================
# Redondear resultados
# ============================================================

columnas_numericas <- sapply(
  resultados_ic,
  is.numeric
)


resultados_ic[
  columnas_numericas
] <- round(
  resultados_ic[
    columnas_numericas
  ],
  digits = 3
)


# ============================================================
# Imprimir tabla final
# ============================================================

print(
  resultados_ic,
  row.names = FALSE
)


# ============================================================
# Tabla reducida para reportar
#
# Se utiliza BCa como intervalo principal.
# ============================================================

tabla_final <- data.frame(
  
  Outcome = c(
    "SBP",
    "DBP"
  ),
  
  ATE = resultados_ic$ATE,
  
  SE_bootstrap = resultados_ic$SE_bootstrap,
  
  IC95_BCa_inferior = resultados_ic$BCa_95_inferior,
  
  IC95_BCa_superior = resultados_ic$BCa_95_superior
)


print(
  tabla_final,
  row.names = FALSE
)


# ============================================================
# Guardar resultados
# ============================================================

write.csv(
  resultados_ic,
  file = file.path(
    ruta_salida,
    "ATE_bootstrap_intervalos.csv"
  ),
  row.names = FALSE
)


write.csv(
  tabla_final,
  file = file.path(
    ruta_salida,
    "ATE_bootstrap_BCa.csv"
  ),
  row.names = FALSE
)





boot_ate_csv <- data.frame(
  iteracion = 1:nrow(boot_ate$t),
  ATE_SBP = boot_ate$t[, 1],
  ATE_DBP = boot_ate$t[, 2]
)

write.csv(
  boot_ate_csv,
  file = file.path(
    ruta_salida,
    "boot_ate_t.csv"
  ),
  row.names = FALSE
)