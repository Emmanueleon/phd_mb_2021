#' @title: Importación datos pools de Qiime
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Importar y fusionar los artefactos generados por QIIME2 (.qza)
#'   de los 10 pools de secuenciación en un objeto phyloseq conjunto (phy_raw).
#'   Incluye metadatos de microbiota, tablas de OTUs, taxonomía y árbol filogenético.


# Entradas:
#   out/qiime2/Pool*/  — artefactos .qza por pool (table, rooted-tree, taxonomy)
#   meta/metadata_Pool*.tsv — metadatos por pool
#   Google Sheets       — metadata conjunto

# Salidas:
#   phy_raw.RDS         — Objeto phyloseq crudo conjunto
#   metadata_mb.RDS     — Metadatos de microbiota
#   otu_singles.RData   — Tablas OTU por pool
#   taxa_singles.RData  — Tablas de taxonomía por pool
#   pools_singles.RData — Objetos phyloseq por pool

# NOTA ÁRBOL FILOGENÉTICO:
#   El árbol actual es provisional (rtree + set.seed).
#   Pendiente sustituir por árbol real de QIIME2 cuando se recuperen
#   los archivos .qza completos con la asesora de tesis.
#   El resto del script NO cambia al hacer esa sustitución.


# Anterior: 03_master_data.R

# Siguiente: 05_calidad_filtrado.R

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
