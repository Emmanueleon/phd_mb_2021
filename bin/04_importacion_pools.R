#' @title: Importación datos pools de Qiime
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Importar y fusionar los artefactos generados por QIIME2 (.biom, .tsv, .nwk) 
#' de los 10 pools de secuenciación en un objeto phyloseq conjunto (phy_raw), sin depender de qiime2R. 



# Entradas:
#   feature-table.biom  — tabla de ASVs
#   taxonomy.tsv        — clasificación taxonómica
#   tree.nwk            — árbol filogenético enraizado
#   meta/metadata_Pool*.tsv — metadatos por pool
#   data/processed/master_data.RDS — metadatos clínicos (opcional, para join)

# Salidas:
#   phy_raw.RDS         — Objeto phyloseq crudo conjunto
#   metadata_mb.RDS     — Metadatos de microbiota
#   pools_singles.RData — Objetos phyloseq por pool


# Anterior: 00_export_qiime2.sh

# Siguiente: 05_calidad_filtrado.R

# ::::::::::::::::::::::::::::
# Preparación del ambiente
# :::::::::::::::::::::::::::

## Limpiar el ambiente
rm(list = ls())

## Cargar librerías
suppressPackageStartupMessages({
    library(phyloseq)
    library(biomformat)
    library(ape)
    library(tidyverse)
    library(here)
})

## Cear carpeta de salida
dir.create(here("out", "phyloseq"), recursive = TRUE, showWarnings = FALSE)

## Rutas de salida
out <- list(
    phy_raw       = here("out", "phyloseq", "phy_raw.RDS"),
    metadata_mb   = here("out", "phyloseq", "metadata_mb.RDS"),
    pools_singles = here("out", "phyloseq", "pools_singles.RData")
)

## Llamado funciones microbiota 
source(here("bin", "functions_microbiota.R"))


# ::::::::::::::::::::::::::::::::::::::::::::::::
# Construir phyloseq por pool
# ::::::::::::::::::::::::::::::::::::::::::::::::

message("Construyendo objetos phyloseq por pool...")

## Primera corrida :::::::::::::::::::::::::::::::::::::::::::::

### Pool A
message("  Pool A...")
phy_a <- build_phyloseq("PoolA", here("meta", "metadata_PoolA.tsv"))
phy_a <- subset_samples(phy_a, sample_id %in% c("B116", "M116", "M120"))
message("    ", nsamples(phy_a), " muestras conservadas")

### Pool B
message("  Pool B...")
phy_b <- build_phyloseq("PoolB", here("meta", "metadata_PoolB.tsv"))
phy_b <- subset_samples(phy_b, sample_id %in% c("B124", "B131", "L134", "M123"))
message("    ", nsamples(phy_b), " muestras conservadas")

### Pool C
message("  Pool C...")
phy_c <- build_phyloseq("PoolC", here("meta", "metadata_PoolC.tsv"))
phy_c <- subset_samples(phy_c, sample_id %in% c("B68", "L140", "L141", "M66"))
message("    ", nsamples(phy_c), " muestras conservadas")

### Pool D
message("  Pool D...")
phy_d <- build_phyloseq("PoolD", here("meta", "metadata_PoolD.tsv"))
phy_d <- subset_samples(phy_d, sample_id %in% c(
    "B118", "B130", "B87",
    "C_Heces", "C_Leche", "H2O",
    "M114", "M122", "M130", "M87", "M94", "MOCK"
))
message("    ", nsamples(phy_d), " muestras conservadas (incluye controles)")
