#' @title: Calidad lecturas
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Analisis de calidad de las lecturas antes, durante y posterior
#'   al proceso de filtrado. Genera tablas comparativas de lecturas crudas vs
#'   filtradas para los tres tipos de muestra (mom, milk, newborn) en un solo
#'   pipeline con group_by() en lugar de scripts separados por muestra.

# Entradas:
#   out/phyloseq/phy_raw.RDS
#   out/phyloseq/phy_filtered.RDS
#   out/phyloseq/phy_objects.RData
#   out/phyloseq/phy_decontam.RData
#   data/processed/master_data.RDS

# Salidas (en out/tables/):
#   reads_summary.RDS   - tabla resumen por tipo de muestra
#   reads_summary.csv   - misma tabla en formato CSV
#   reads_detail.RDS    - tabla detallada por muestra individual
#   reads_detail.csv    - misma tabla en formato CSV

# Anterior: 05_calidad_filtrado.R
# Siguiente: 07_composicion_core.R

# ::::::::::::::::::::::::::::
# Preparacion del ambiente
# :::::::::::::::::::::::::::
rm(list = ls())

suppressPackageStartupMessages({
    library(phyloseq)
    library(tidyverse)
    library(here)
})

## Crear directorio de salida
dir.create(here("out", "tables"), recursive = TRUE, showWarnings = FALSE)

## Rutas de salida
out <- list(
    reads_summary = here("out", "tables", "reads_summary.RDS"),
    reads_summary_csv = here("out", "tables", "reads_summary.csv"),
    reads_detail  = here("out", "tables", "reads_detail.RDS"),
    reads_detail_csv  = here("out", "tables", "reads_detail.csv")
)

## Cargar funciones
source(here("bin", "funciones_microbiota.R"))

## Cargar datos
message("Cargando datos...")
phy_raw      <- readRDS(here("out", "phyloseq", "phy_raw.RDS"))
phy_filtered <- readRDS(here("out", "phyloseq", "phy_filtered.RDS"))
load(here("out", "phyloseq", "phy_objects.RData"))
load(here("out", "phyloseq", "phy_decontam.RData"))

## Variables
sample_types <- c("mom", "milk", "newborn")
taxa_rank    <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1. Resumen de secuencias crudas
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Secuencias crudas ---")

## Distribucion por tipo de muestra
message("Distribucion de muestras:")
print(table(sample_data(phy_raw)$sample))

## Total de lecturas
reads_raw <- reads_tbl(phy_raw, sample_levels = c("mom", "milk", "newborn",
                                            "positive", "negative"))
message("Total lecturas crudas: ", sum(reads_raw$reads))

## Riqueza taxonomica cruda
message("\nRiqueza taxonomica cruda:")
for (tax in taxa_rank) {
    n <- tax_table(phy_raw)[, tax] %>% unique() %>% na.exclude() %>% length()
    message("  ", tax, ": ", n)
}


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Resumen del proceso de filtrado
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Resumen del filtrado ---")

message("Decontam lote 1 — ASVs removidos: ",
        ntaxa(phy_first) - ntaxa(phy_decontam_first))
message("Decontam lote 2 — ASVs removidos: ",
        ntaxa(phy_second) - ntaxa(phy_decontam_second))
message("Muestras < 5000 reads — removidas: ",
        nsamples(phy_raw) - nsamples(phy_filsam))
message("ASVs sin lecturas — removidos: ",
        ntaxa(phy_filsam) - ntaxa(phy_filreadtax))
message("Cloroplastos/Mitocondria — removidos: ",
        ntaxa(phy_filreadtax) - ntaxa(phy_filtax))
message("Controles — removidos: ",
        nsamples(phy_filtax) - nsamples(phy_filtered))


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3. Tabla detallada por muestra: crudas vs filtradas
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Tabla detallada por muestra ---")

## Lecturas crudas
reads_raw_tbl <- reads_tbl(
    phy_raw,
    sample_levels = c("mom", "milk", "newborn", "positive", "negative")
) %>%
    rename(raw_reads = reads)

## Lecturas filtradas
reads_filtered_tbl <- reads_tbl(
    phy_filtered,
    sample_levels = sample_types
) %>%
    rename(filtered_reads = reads) %>%
    select(sample_id, filtered_reads)

## Tabla comparativa
reads_detail <- reads_raw_tbl %>%
    full_join(reads_filtered_tbl, by = "sample_id") %>%
    mutate(
        filtered_reads = replace_na(filtered_reads, 0),
        loss_reads     = raw_reads - filtered_reads,
        pct_retained   = round(filtered_reads / raw_reads * 100, 1),
        status         = case_when(
            filtered_reads == 0 ~ "Excluida",
            TRUE                ~ "Incluida"
        )
    ) %>%
    arrange(sample, desc(raw_reads))

## Muestras excluidas
excluded <- reads_detail %>% filter(status == "Excluida")
message("Muestras excluidas (", nrow(excluded), "):")
print(excluded %>% select(sample_id, sample, raw_reads))


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4. Resumen estadistico por tipo de muestra
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Resumen por tipo de muestra ---")

reads_summary <- reads_detail %>%
    filter(sample %in% sample_types) %>%
    group_by(sample) %>%
    summarize(
        n_total        = n(),
        n_incluidas    = sum(status == "Incluida"),
        n_excluidas    = sum(status == "Excluida"),
        sum_raw        = sum(raw_reads),
        sum_filtered   = sum(filtered_reads),
        mean_raw       = round(mean(raw_reads), 0),
        mean_filtered  = round(mean(filtered_reads[filtered_reads > 0]), 0),
        median_filtered = round(median(filtered_reads[filtered_reads > 0]), 0),
        min_filtered   = min(filtered_reads[filtered_reads > 0]),
        max_filtered   = max(filtered_reads[filtered_reads > 0]),
        .groups = "drop"
    )

print(reads_summary)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 5. Riqueza taxonomica post-filtrado
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Riqueza taxonomica post-filtrado ---")

message("Distribucion de muestras filtradas:")
print(table(sample_data(phy_filtered)$sample))

reads_filtered_total <- reads_tbl(phy_filtered, sample_levels = sample_types)
message("Total lecturas filtradas: ", sum(reads_filtered_total$reads))

for (tax in taxa_rank) {
    n <- tax_table(phy_filtered)[, tax] %>% unique() %>% na.exclude() %>% length()
    message("  ", tax, ": ", n)
}


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 6. Exportacion
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Exportacion ---")

saveRDS(reads_summary, out$reads_summary)
write.csv(reads_summary, out$reads_summary_csv, row.names = FALSE)
message("reads_summary guardado")

saveRDS(reads_detail, out$reads_detail)
write.csv(reads_detail, out$reads_detail_csv, row.names = FALSE)
message("reads_detail guardado")

message("\n--- 06_calidad_lecturas.R completado ---")
message("Siguiente paso: 07_composicion_core.R")