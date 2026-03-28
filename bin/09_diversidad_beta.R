#' @title: Diversidad beta
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Analisis de diversidad beta usando distancias UniFrac no
#'   ponderada (riqueza) y ponderada (abundancia) con PCoA y pruebas
#'   estadisticas (PERMANOVA, Betadisper, ANOSIM) para los tres tipos de
#'   muestra y cinco variables clinicas.

# Entradas:
#   out/phyloseq/phy_filtered.RDS

# Salidas:
#   beta_{sample}_{variable}_{distancia}.RDS — resultados por combinacion

# Anterior: 08_diversidad_alfa.R
# Siguiente: 10_lefse.R

# ::::::::::::::::::::::::::::
# Preparacion del ambiente
# :::::::::::::::::::::::::::
rm(list = ls())

suppressPackageStartupMessages({
    library(phyloseq)
    library(microbiome)
    library(vegan)
    library(tidyverse)
    library(here)
})

## Crear directorio de salida
dir.create(here("out", "tables", "beta"), recursive = TRUE, showWarnings = FALSE)

## Cargar funciones
source(here("bin", "funciones_microbiota.R"))

## Cargar datos
message("Cargando datos...")
phy_filtered <- readRDS(here("out", "phyloseq", "phy_filtered.RDS"))
message("phy_filtered: ", nsamples(phy_filtered), " muestras")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Analisis del articulo
# Combinaciones reportadas en Cervantes-Monroy et al. 2024
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Analisis del articulo ---")

## IMC y adiposidad materna
beta_mom_imc_u   <- run_beta(phy_filtered, "mom",     "imc",      "unifrac")
beta_mom_imc_w   <- run_beta(phy_filtered, "mom",     "imc",      "wunifrac")
beta_mom_fat_u   <- run_beta(phy_filtered, "mom",     "fat",      "unifrac")
beta_mom_fat_w   <- run_beta(phy_filtered, "mom",     "fat",      "wunifrac")

## Edad materna
beta_mom_age_u   <- run_beta(phy_filtered, "mom",     "mom_age30", "unifrac")
beta_mom_age_w   <- run_beta(phy_filtered, "mom",     "mom_age30", "wunifrac")

## Modo de parto
beta_mom_del_u   <- run_beta(phy_filtered, "mom",     "delivery", "unifrac")
beta_mom_del_w   <- run_beta(phy_filtered, "mom",     "delivery", "wunifrac")
beta_milk_del_u  <- run_beta(phy_filtered, "milk",    "delivery", "unifrac")
beta_milk_del_w  <- run_beta(phy_filtered, "milk",    "delivery", "wunifrac")
beta_nb_del_u    <- run_beta(phy_filtered, "newborn", "delivery", "unifrac")
beta_nb_del_w    <- run_beta(phy_filtered, "newborn", "delivery", "wunifrac")

## Sexo del recien nacido
beta_mom_sex_u   <- run_beta(phy_filtered, "mom",     "sex",      "unifrac")
beta_mom_sex_w   <- run_beta(phy_filtered, "mom",     "sex",      "wunifrac")
beta_milk_sex_u  <- run_beta(phy_filtered, "milk",    "sex",      "unifrac")
beta_milk_sex_w  <- run_beta(phy_filtered, "milk",    "sex",      "wunifrac")
beta_nb_sex_u    <- run_beta(phy_filtered, "newborn", "sex",      "unifrac")
beta_nb_sex_w    <- run_beta(phy_filtered, "newborn", "sex",      "wunifrac")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Exportacion
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Exportacion ---")

out_dir <- here("out", "tables", "beta")

resultados <- list(
    beta_mom_imc_u,  beta_mom_imc_w,
    beta_mom_fat_u,  beta_mom_fat_w,
    beta_mom_age_u,  beta_mom_age_w,
    beta_mom_del_u,  beta_mom_del_w,
    beta_milk_del_u, beta_milk_del_w,
    beta_nb_del_u,   beta_nb_del_w,
    beta_mom_sex_u,  beta_mom_sex_w,
    beta_milk_sex_u, beta_milk_sex_w,
    beta_nb_sex_u,   beta_nb_sex_w
)

walk(resultados, function(r) {
    nombre <- paste0("beta_", r$sample_type, "_",
                     r$variable, "_", r$distancia)
    saveRDS(r, file.path(out_dir, paste0(nombre, ".RDS")))
    message("  ", nombre, " guardado")
})

message("\n--- 09_diversidad_beta.R completado ---")
message("Archivos en: ", out_dir)
message("Siguiente paso: 10_lefse.R")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# NOTA: Como anadir nuevas combinaciones
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Para anadir una combinacion nueva, simplemente llama run_beta()
# con los argumentos deseados y anadela a la lista de resultados:
#
#   # Ejemplo: sexo del bebe en microbiota materna
#   beta_mom_sex_u <- run_beta(phy_filtered, "mom", "sex", "unifrac")
#
#   # Ejemplo: edad materna en leche
#   beta_mom_age_u <- run_beta(phy_filtered, "mom", "mom_age30", "unifrac")
#
# Para correr TODAS las combinaciones posibles de una vez:
#
#   sample_types <- c("mom", "milk", "newborn")
#   variables    <- c("imc", "fat", "delivery", "sex", "mom_age30")
#   distancias   <- c("unifrac", "wunifrac")
#
#   combinaciones <- expand.grid(
#       sample_type = sample_types,
#       variable    = variables,
#       distancia   = distancias,
#       stringsAsFactors = FALSE
#   )
#
#   todos_resultados <- pmap(combinaciones, function(sample_type, variable, distancia) {
#       run_beta(phy_filtered, sample_type, variable, distancia)
#   })
#
# ADVERTENCIA: expand.grid genera 30 combinaciones (3x5x2).
# Algunas pueden tener n insuficiente — revisar warnings antes de reportar.
