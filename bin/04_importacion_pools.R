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

