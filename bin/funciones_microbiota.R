#' @title: Funciones para el análisis de microbiota
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Colección de funciones reutilizables para el análisis de
#'   microbiota. Se carga al inicio de cada script de análisis con:
#'   source(here("bin", "functions_microbiota.R"))

# Paquetes requeridos
# library(tidyverse), library(phyloseq), library(microbiome),
# library(biomformat), library(ape), library(broom), library(Hmisc)

# =============================================================================
# 1. IMPORTACIÓN PHYLOSEQ
# =============================================================================

#' build_phyloseq
#' @description Construye un objeto phyloseq desde los archivos exportados por QIIME2
#' @param pool Nombre del pool (ej. "PoolA")
#' @param metadata_path Ruta al archivo .tsv de metadatos del pool
#' @return Objeto phyloseq con OTU, taxonomía, árbol y metadatos

build_phyloseq <- function(pool, metadata_path) {
    export_dir <- here::here("out", "export", pool)

    ## Feature table
    biom_data <- biomformat::read_biom(file.path(export_dir, "feature-table.biom"))
    otu_matrix <- biomformat::biom_data(biom_data) %>% as.matrix()
    otu <- phyloseq::otu_table(otu_matrix, taxa_are_rows = TRUE)

    ## Taxonomía
    tax_raw <- readr::read_tsv(
        file.path(export_dir, "taxonomy.tsv"),
        show_col_types = FALSE
    )

    tax_mat <- tax_raw %>%
        tidyr::separate(Taxon,
            into = c(
                "Kingdom", "Phylum", "Class", "Order",
                "Family", "Genus", "Species"
            ),
            sep = ";\\s*", fill = "right"
        ) %>%
        tibble::column_to_rownames("Feature ID") %>%
        dplyr::select(-Confidence) %>%
        as.matrix()

    tax <- phyloseq::tax_table(tax_mat)
    tree <- ape::read.tree(file.path(export_dir, "tree.nwk"))

    meta_raw <- readr::read_tsv(
        metadata_path,
        show_col_types = FALSE,
        comment = "#"
    ) %>%
        dplyr::select(-`sample-id`) %>% # eliminar columna duplicada
        tibble::column_to_rownames("sample_id")

    meta <- phyloseq::sample_data(meta_raw)
    phyloseq::phyloseq(otu, tax, tree, meta)
}

message("funciones_microbiota.R cargado — build_phyloseq disponible")
