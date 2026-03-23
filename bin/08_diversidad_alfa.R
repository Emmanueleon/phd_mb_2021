#' @title: Diversidad alfa
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Analisis de diversidad alfa usando metricas de riqueza,
#'   diversidad, equitatividad y Faith PD para los tres tipos de muestra.
#'   Se usa el objeto rarefaccionado a 5000 lecturas.

# Entradas:
#   out/phyloseq/phy_objects.RData  (contiene phy_rarefied5k)

# Salidas (en out/tables/alfa/):
#   adiv_table.RDS      - tabla completa de metricas alfa
#   adiv_table.csv
#   adiv_{sample}.RDS   - tabla por tipo de muestra
#   adiv_{sample}.csv

# Anterior: 07_composicion_core.R
# Siguiente: 09_diversidad_beta.R

# ::::::::::::::::::::::::::::
# Preparacion del ambiente
# :::::::::::::::::::::::::::
rm(list = ls())

suppressPackageStartupMessages({
    library(phyloseq)
    library(microbiome)
    library(picante)
    library(tidyverse)
    library(here)
})

## Crear directorio de salida
dir.create(here("out", "tables", "alfa"), recursive = TRUE, showWarnings = FALSE)

## Rutas de salida
out_dir <- here("out", "tables", "alfa")

## Cargar funciones
source(here("bin", "funciones_microbiota.R"))

## Cargar datos
message("Cargando datos...")
load(here("out", "phyloseq", "phy_objects.RData"))

message("phy_rarefied5k: ", nsamples(phy_rarefied5k), " muestras")
message("Arbol: ", ntaxa(phy_rarefied5k), " tips")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1. Metricas de diversidad alfa
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Calculando metricas alfa ---")

alfa_metrics <- c(
    "Shannon", "Simpson", "InvSimpson", "Fisher",
    "Observed", "Chao1", "ACE"
)

## Tabla completa con todas las muestras
adiv_tbl <- adiv_table(phy_rarefied5k, metric = alfa_metrics)

message(
    "Tabla adiv: ", nrow(adiv_tbl), " muestras, ",
    ncol(adiv_tbl), " variables"
)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Distribucion de metricas por tipo de muestra
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Resumen por tipo de muestra ---")

adiv_summary <- adiv_tbl %>%
    group_by(sample) %>%
    summarize(
        n = n(),
        Shannon_mean = round(mean(Shannon, na.rm = TRUE), 2),
        Shannon_sd = round(sd(Shannon, na.rm = TRUE), 2),
        Chao1_mean = round(mean(Chao1, na.rm = TRUE), 0),
        Chao1_sd = round(sd(Chao1, na.rm = TRUE), 0),
        PD_mean = round(mean(PD, na.rm = TRUE), 2),
        PD_sd = round(sd(PD, na.rm = TRUE), 2),
        .groups = "drop"
    )

print(adiv_summary)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3. Tablas por tipo de muestra
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Generando tablas por tipo de muestra ---")

sample_types <- c("mom", "milk", "newborn")

adiv_by_sample <- map(sample_types, function(s) {
    adiv_tbl %>%
        filter(sample == s) %>%
        select(
            key, sample, no_reads,
            all_of(alfa_metrics), PD,
            pielou, starts_with("fisher")
        )
}) %>% set_names(sample_types)

## Resumen por tipo
walk(sample_types, function(s) {
    message("  ", s, ": ", nrow(adiv_by_sample[[s]]), " muestras")
})


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4. Exportacion
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Exportacion ---")

## Tabla completa
saveRDS(adiv_tbl, file.path(out_dir, "adiv_table.RDS"))
write.csv(adiv_tbl, file.path(out_dir, "adiv_table.csv"), row.names = FALSE)
message("adiv_table guardada")

## Tablas por tipo de muestra
walk(sample_types, function(s) {
    saveRDS(
        adiv_by_sample[[s]],
        file.path(out_dir, paste0("adiv_", s, ".RDS"))
    )
    write.csv(adiv_by_sample[[s]],
        file.path(out_dir, paste0("adiv_", s, ".csv")),
        row.names = FALSE
    )
    message("  adiv_", s, " exportado")
})

message("\n--- 08_diversidad_alfa.R completado ---")
message("Archivos en: ", out_dir)
message("Siguiente paso: 09_diversidad_beta.R")
