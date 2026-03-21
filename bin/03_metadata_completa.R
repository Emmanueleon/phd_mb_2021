#' @title: Construcción Metadatos
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Construir la mega base de datos uniendo todas las tablas individuales
#   (crudas y derivadas) mediante left_join explícito por key
#   Esta base es el punto de partida para todos los análisis estadísticos
#   posteriores. No se modifica aquí — solo se ensambla y valida.


# Entradas:
#   data/raw/      — bases crudas de 01_importacion.R
#   data/derived/  — variables clínicas de 02_variables_clinicas.R

# Salidas:
#   master_data.RDS        — Mega base completa
#   colnames.RDA        — Vectores de nombres de columna por tabla

# Anterior: 02_variables_clinicas.R

# Siguiente: cualquier script de análisis (04 en adelante)

# ::::::::::::::::::::::::::::
# Preparación del ambiente
# :::::::::::::::::::::::::::

## Limpiar el ambiente
rm(list = ls())

## Cargar librerías
suppressPackageStartupMessages({
    library(tidyverse)
    library(here)
})

## Cear carpeta de salida
dir.create(here("data", "processed"), recursive = TRUE, showWarnings = FALSE)


# ::::::::::::::::::::::::::::::
# Cargar todas las bases
# ::::::::::::::::::::::::::::::

## Datos crudos
mom_demo <- readRDS(here("data", "raw", "mom_demo.RDS"))
hm_char <- readRDS(here("data", "raw", "hm_char.RDS"))
hm_practice <- readRDS(here("data", "raw", "hm_practice.RDS"))
feces_char <- readRDS(here("data", "raw", "feces_char.RDS"))

## Datos con ingeniería de variables
mom_anthro <- readRDS(here("data", "derived", "mom_anthro_d.RDS"))
mom_bx <- readRDS(here("data", "derived", "mom_bx_d.RDS"))
mom_diet <- readRDS(here("data", "derived", "mom_diet_d.RDS"))
nb_anthro <- readRDS(here("data", "derived", "nb_anthro_d.RDS"))

## Registro del n de partida
n_inicio <- nrow(mom_demo)
message("Participantes en mom_demo (base de referencia): ", n_inicio)


# ::::::::::::::::::::::::::::::
# Construcción de la mega base
# ::::::::::::::::::::::::::::::
master_data <- mom_demo %>%
    left_join(mom_anthro, by = "key") %>%
    left_join(mom_bx, by = "key") %>%
    left_join(mom_diet, by = "key") %>%
    left_join(nb_anthro, by = "key") %>%
    left_join(hm_char, by = "key") %>%
    left_join(hm_practice, by = "key") %>%
    left_join(feces_char, by = "key")

# -----------------------------------------------------------------------------
# Validación
# -----------------------------------------------------------------------------
n_final <- nrow(master_data)

message("\n--- Validación del join ---")
message("Filas antes : ", n_inicio)
message("Filas después: ", n_final)

if (n_final != n_inicio) {
    warning("El número de filas cambió tras el join. Revisar duplicados en px_num.")
} else {
    message("OK — n se conservó correctamente.")
}

## Revisar duplicados en la llave
duplicados <- master_data %>%
    count(key) %>%
    filter(n > 1)
if (nrow(duplicados) > 0) {
    warning("Se detectaron px_num duplicados: ")
    print(duplicados)
} else {
    message("OK — px_num único por participante.")
}

## Reporte de completitud por tabla
message("\n--- Completitud de variables clave ---")
vars_clave <- c(
    "age.calculated", "imc.2", "bd.fat.2", # mom_anthro
    "delivery.mode", "gender", # nb_anthro
    "iom.gain", "mom_age30", "nutc.imc.2", # derivadas
    "glucose", "cholesterol", # mom_bx
    "hm.1"
) # hm_char
for (v in vars_clave) {
    if (v %in% colnames(master_data)) {
        n_na <- sum(is.na(master_data[[v]]))
        message("  ", v, ": ", n_na, " NAs (", round(n_na / n_final * 100, 1), "%)")
    }
}

# ::::::::::::::::::::::::::::::
# Exportación
# ::::::::::::::::::::::::::::::
saveRDS(master_data, here("data", "processed", "master_data.RDS"))
message("\nmaster_data guardado: ", nrow(master_data), " filas x ", ncol(master_data), " columnas")

## Vectores de nombres de columna por tabla (útiles para selección posterior)
colnames_mom_demo <- colnames(mom_demo)
colnames_mom_anthro <- colnames(mom_anthro)
colnames_mom_bx <- colnames(mom_bx)
colnames_mom_diet <- colnames(mom_diet)
colnames_nb_anthro <- colnames(nb_anthro)
colnames_hm_char <- colnames(hm_char)
colnames_hm_practice <- colnames(hm_practice)
colnames_feces_char <- colnames(feces_char)

save(colnames_mom_demo,
    colnames_mom_anthro,
    colnames_mom_bx,
    colnames_mom_diet,
    colnames_nb_anthro,
    colnames_hm_char,
    colnames_hm_practice,
    colnames_feces_char,
    file = here("data","processed","colnames.RDA")
)

message("colnames.RDA guardado.")

# -----------------------------------------------------------------------------
# Ejemplo de uso para análisis posteriores
# -----------------------------------------------------------------------------
# Para seleccionar solo las variables necesarias en un análisis específico:
#
#   load(here("data", "colnames.RDA"))
#   master_data <- readRDS(here("data", "master_data.RDS"))
#
#   # Ejemplo: IMC materno vs antropometría neonatal
#   datos_analisis <- master_data %>%
#     select(px_num, nutc.imc.2, nutc.fat.2,
#            all_of(colnames_nb_anthro),
#            z.birth.wa, z.nb.wa)
#
# master_data nunca se modifica — cada análisis crea su propio subconjunto.

message("\n--- 03_mega_base.R completado ---")
message("Siguiente paso: scripts de análisis (04 en adelante)")
