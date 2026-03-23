#' @title: Análisis composicional y core microbiota de la triada
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Analisis de la composicion taxonomica y core microbiota para
#'   los tres tipos de muestra (mom, milk, newborn) en un solo pipeline.
#'   Genera tablas de abundancia relativa con etiquetado de taxones prevalentes
#'   y core microbiota al 70% de prevalencia.


# Entradas:
# out/phyloseq/phy_filtered.RDS

# Salidas:
#   metataxa_prevalence_{sample}.RDS  - tabla larga con prevalencia por taxon
#   abundance_table_{sample}.RDS      - tabla de abundancia relativa
#   abundance_table_{sample}.csv
#   abundance_table_final_{sample}.csv - con filo, orden, familia
#   tax_table_{sample}.RDS            - tabla taxonomica por muestra

# Anterior: 06_calidad_lecturas.R
# Siguiente: 08_diversidad_alfa.R

# ::::::::::::::::::::::::::::
# Preparacion del ambiente
# :::::::::::::::::::::::::::
rm(list = ls())

suppressPackageStartupMessages({
    library(phyloseq)
    library(microbiome)
    library(tidyverse)
    library(here)
})

## Crear directorio de salida
dir.create(here("out", "tables", "composicion"), recursive = TRUE, showWarnings = FALSE)

## Cargar funciones
source(here("bin", "funciones_microbiota.R"))

## Cargar datos
message("Cargando datos...")
phy_filtered <- readRDS(here("out", "phyloseq", "phy_filtered.RDS"))

 
message("phy_filtered: ", nsamples(phy_filtered), " muestras, ",
        ntaxa(phy_filtered), " ASVs")
 
## Taxones no descriptivos a excluir de las listas de prevalencia
uncertain_taxa <- c("Unclassified", "Incertae_Sedis", "uncultured")
 
## Tipos de muestra a analizar
sample_types <- c("mom", "milk", "newborn")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Funcion principal por tipo de muestra
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
analyze_composition <- function(sample_type_val, phy_filtered, master_data) {
 
    message("\n--- Procesando: ", sample_type_val, " ---")
 
    # ── 1. Objeto composicional por tipo de muestra ──────────────
    phy_comp <- prune_samples(
    sample_data(phy_filtered)$sample == sample_type_val,
    phy_filtered) %>%
        microbiome::transform(transform = "compositional")
 
    message("  Muestras: ", nsamples(phy_comp))
 
    # ── 2. Tabla taxonomica larga ─────────────────────────────────
    taxa_tbl <- taxa_table(phy_comp)
 
    # ── 3. Union con metadatos clinicos ──────────────────────────
    metadata_taxa <- taxa_tbl
 
    # ── 4. Listas de taxones prevalentes ─────────────────────────
    ### Filos top 4
    phylum_lst <- metadata_taxa %>%
        top_taxa_lst(tax_rank = Phylum, top = 4)
 
    message("  Filos top 4: ", paste(phylum_lst, collapse = ", "))
 
    ### Familias top 5 por filo (excluyendo no descriptivos)
    family_lst <- metadata_taxa %>%
        filter(!Family %in% uncertain_taxa) %>%
        taxa_lst(tax_rank = Family, top = 5, phylum_vector = phylum_lst)
 
    ### Generos top 5 por filo (excluyendo no descriptivos)
    genus_lst <- metadata_taxa %>%
        filter(!Genus %in% uncertain_taxa) %>%
        taxa_lst(tax_rank = Genus, top = 5, phylum_vector = phylum_lst)
 
    # ── 5. Etiquetado de taxones prevalentes ─────────────────────
    metadata_taxa_prevalence <- metadata_taxa %>%
        mutate(
            prevalence_phylum = case_when(
                Phylum %in% phylum_lst ~ Phylum,
                TRUE ~ "Other"
            ),
            prevalence_family = case_when(
                Family %in% family_lst ~ Family,
                TRUE ~ "Other"
            ),
            prevalence_genus = case_when(
                Genus %in% genus_lst ~ Genus,
                TRUE ~ "Other"
            )
        ) %>%
        ## Etiquetas "Other X" por filo para familia
        mutate(prevalence_family = case_when(
            prevalence_phylum == "Actinobacteriota" & prevalence_family == "Other" ~ "Other Actinobacteriota",
            prevalence_phylum == "Firmicutes"       & prevalence_family == "Other" ~ "Other Firmicutes",
            prevalence_phylum == "Proteobacteria"   & prevalence_family == "Other" ~ "Other Proteobacteria",
            prevalence_phylum == "Bacteroidota"     & prevalence_family == "Other" ~ "Other Bacteroidota",
            prevalence_phylum == "Other"                                            ~ "Other",
            TRUE ~ Family
        )) %>%
        ## Etiquetas "Other X" por filo para genero
        mutate(prevalence_genus = case_when(
            prevalence_phylum == "Actinobacteriota" & prevalence_genus == "Other" ~ "Other Actinobacteriota",
            prevalence_phylum == "Firmicutes"       & prevalence_genus == "Other" ~ "Other Firmicutes",
            prevalence_phylum == "Proteobacteria"   & prevalence_genus == "Other" ~ "Other Proteobacteria",
            prevalence_phylum == "Bacteroidota"     & prevalence_genus == "Other" ~ "Other Bacteroidota",
            prevalence_phylum == "Other"                                           ~ "Other",
            TRUE ~ Genus
        )) %>%
        mutate(across(Kingdom:prevalence_genus, as.factor)) %>%
        mutate(
            prevalence_phylum = factor(prevalence_phylum,
                                       levels = c(phylum_lst, "Other")),
            prevalence_family = fct_reorder(prevalence_family,
                                            10 * as.integer(prevalence_phylum) +
                                            grepl("Other", prevalence_family)),
            prevalence_genus  = fct_reorder(prevalence_genus,
                                            10 * as.integer(prevalence_phylum) +
                                            grepl("Other", prevalence_genus))
        )
 
    # ── 6. Tabla de abundancia relativa ──────────────────────────
    abundance_tbl <- taxa_relabun_by_group(
        df        = metadata_taxa_prevalence,
        var       = sample,
        taxa_fill = prevalence_phylum,
        group_var = sample_type_val,
        taxa_rank = Genus
    )$by_group
 
    message("  Generos en tabla de abundancia: ", nrow(abundance_tbl))
    message("  Suma total (debe ser ~100): ", round(sum(abundance_tbl$suma), 1))
 
    # ── 7. Core microbiota al 70% ─────────────────────────────────
    phy_sample <- prune_samples(
    sample_data(phy_filtered)$sample == sample_type_val,
    phy_filtered)
    core_tbl   <- core_taxa(phy_sample, preval = 0.70, tax_rank = Genus)
    core_lst   <- pull(core_tbl, Genus)
 
    message("  Generos core (70%): ", length(core_lst))
    message("  ", paste(core_lst, collapse = ", "))
 
    ## Añadir columna core a la tabla de abundancia
    abundance_tbl <- abundance_tbl %>%
        mutate(core = case_when(
            Genus %in% core_lst ~ "Yes",
            TRUE ~ "No"
        ))
 
    # ── 9. Tabla final con jerarquia taxonomica ───────────────────
    taxa_cols <- taxa_tbl %>%
        select(Phylum, Order:Genus) %>%
        distinct(Genus, .keep_all = TRUE)
 
    abundance_final <- abundance_tbl %>%
        inner_join(taxa_cols, by = "Genus") %>%
        select(sample, Phylum:Family, Genus:core)
 
    # ── 10. Generos > 1% abundancia ───────────────────────────────
    genus_1pct <- abundance_tbl %>%
        filter(suma >= 1) %>%
        mutate(Genus = as.character(Genus)) %>%
        pull(Genus)
 
    message("  Generos > 1% abundancia: ", length(genus_1pct))
 
    ## Retornar todos los objetos
    list(
        metadata_taxa_prevalence = metadata_taxa_prevalence,
        abundance_tbl            = abundance_tbl,
        abundance_final          = abundance_final,
        core_tbl                 = core_tbl,
        core_lst                 = core_lst,
        genus_1pct               = genus_1pct,
        phylum_lst               = phylum_lst,
        family_lst               = family_lst,
        genus_lst                = genus_lst
    )
}
 
 
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Correr para los 3 tipos de muestra
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
results <- map(sample_types, ~ analyze_composition(
    sample_type_val = .x,
    phy_filtered    = phy_filtered,
    master_data     = master_data
)) %>% set_names(sample_types)
 
 
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Exportacion
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Exportacion ---")
 
out_dir <- here("out", "tables", "composicion")
 
walk(sample_types, function(s) {
 
    r <- results[[s]]
 
    saveRDS(r$metadata_taxa_prevalence,
            file.path(out_dir, paste0("metataxa_prevalence_", s, ".RDS")))
 
    saveRDS(r$abundance_tbl,
            file.path(out_dir, paste0("abundance_table_", s, ".RDS")))
 
    write.csv(r$abundance_tbl,
              file.path(out_dir, paste0("abundance_table_", s, ".csv")),
              row.names = FALSE)
 
    write.csv(r$abundance_final,
              file.path(out_dir, paste0("abundance_table_final_", s, ".csv")),
              row.names = FALSE)

 
    message("  ", s, " exportado")
})
 
message("\n--- 07_composicion_core.R completado ---")
message("Archivos en: ", out_dir)
message("Siguiente paso: 08_diversidad_alfa.R")
 

