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

#' @title build_phyloseq
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


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Calidad y filtrado
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#' @title reads_tbl
#' @description Crear una tabla con el número de lecturas de secuenciación por muestra para determinado objeto phyloseq.
#' @param phy_object Objeto phyloseq a evaluar
#' @param sample_levels Vector de niveles de muestra a considerar
#' @return Tabla con el número de lecturas por muestra

reads_tbl <- function(phy_object, sample_levels) {
    sample_data(phy_object)$reads <- sample_sums(phy_object)
    data.frame(sample_data(phy_object)) %>%
        tibble::rownames_to_column(var = "sample_id") %>% # convertir rownames a columna para formato tidyverse
        as_tibble() %>%
        filter(sample %in% sample_levels) %>%
        select(sample_id, sample, reads)
}

# :::::::::::::::::::::::::::::::::::::::::
# Composicion y core microbiota
# :::::::::::::::::::::::::::::::::::::::::
#' @title taxa_table
#' @description Crea tabla larga con abundancias y rangos taxonomicos por muestra
#' @param phy_object Objeto phyloseq (composicional)
#' @return Tibble con sample_id, Abundance, OTU, Kingdom:Species
taxa_table <- function(phy_object) {
    phy_object %>%
        phyloseq::psmelt() %>%
        tibble::as_tibble() %>%
        dplyr::rename(key = Sample) %>%
        dplyr::select(key, sample, Abundance, OTU, Kingdom:Species) %>%
        dplyr::mutate(across(Kingdom:Species,
                             ~ gsub("^[a-z]__", "", .x))) %>%
        tidyr::replace_na(list(
            Kingdom = "Unclassified", Phylum  = "Unclassified",
            Class   = "Unclassified", Order   = "Unclassified",
            Family  = "Unclassified", Genus   = "Unclassified",
            Species = "Unclassified"
        ))
}

#' @title top_taxa_lst
#' @description Identifica los top n taxones mas abundantes
#' @param df_taxa Tabla de taxones (output de taxa_table)
#' @param tax_rank Rango taxonomico (unquoted, ej. Genus)
#' @param top Numero de taxones a identificar
#' @return Vector con los nombres de los top taxones
top_taxa_lst <- function(df_taxa, tax_rank, top) {
    df_taxa %>%
        dplyr::select({{ tax_rank }}, Abundance) %>%
        dplyr::group_by({{ tax_rank }}) %>%
        dplyr::summarise(sum_abun = sum(Abundance), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(sum_abun)) %>%
        dplyr::distinct({{ tax_rank }}) %>%
        dplyr::slice_head(n = top) %>%
        dplyr::pull()
}

#' @title taxa_lst
#' @description Identifica los top n taxones por filo
#' @param df_taxa Tabla de taxones
#' @param tax_rank Rango taxonomico (unquoted)
#' @param top Numero de taxones por filo
#' @param phylum_vector Vector de filos de interes
#' @return Vector con nombres de taxones
taxa_lst <- function(df_taxa, tax_rank, top, phylum_vector) {
    df_taxa %>%
        dplyr::select(Phylum, {{ tax_rank }}, Abundance) %>%
        dplyr::group_by(Phylum, {{ tax_rank }}) %>%
        dplyr::summarise(sum_abundance = sum(Abundance), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(sum_abundance)) %>%
        dplyr::group_by(Phylum) %>%
        dplyr::slice(1:top) %>%
        dplyr::ungroup() %>%
        dplyr::filter(Phylum %in% phylum_vector) %>%
        dplyr::pull({{ tax_rank }})
}

#' @title taxa_relabun_by_group
#' @description Tabla de abundancia relativa promedio por grupo y rango taxonomico
#' @param df Tabla larga con taxones
#' @param var Variable de agrupacion (unquoted)
#' @param taxa_fill Rango taxonomico para colorear (unquoted)
#' @param group_var Valor del grupo a filtrar (string)
#' @param taxa_rank Rango taxonomico para agrupar (unquoted)
#' @return Lista con dos tablas: abundances (por muestra) y by_group (por grupo)
taxa_relabun_by_group <- function(df, var, taxa_fill, group_var, taxa_rank) {
    samtaxa_sum <- df %>%
        dplyr::group_by({{ var }}, key, Phylum, Family, Genus, {{ taxa_fill }}) %>%
        dplyr::summarise(rel_abund = sum(Abundance), .groups = "drop") %>%
        dplyr::group_by({{ var }}, Phylum, Family, Genus, {{ taxa_fill }}) %>%
        dplyr::summarise(mean_rel_abun = 100 * mean(rel_abund), .groups = "drop")

    var_abundance <- samtaxa_sum %>%
        dplyr::group_by({{ var }}, {{ taxa_rank }}) %>%
        dplyr::summarise(suma = sum(mean_rel_abun), .groups = "drop") %>%
        dplyr::arrange({{ var }}, dplyr::desc(suma)) %>%
        dplyr::filter(suma != 0 & {{ var }} == group_var)

    list(abundances = samtaxa_sum, by_group = var_abundance)
}


#' @title core_taxa
#' @description Tabla de abundancia de taxones del core microbiota
#' @param phy_object Objeto phyloseq
#' @param preval Umbral de prevalencia (ej. 0.70 para 70 por ciento)
#' @param tax_rank Rango taxonomico (unquoted)
#' @return Tibble con taxones del core y su abundancia total
core_taxa <- function(phy_object, preval, tax_rank) {
    core_list   <- microbiome::core_members(phy_object, prevalence = preval)
    sample_core <- phyloseq::prune_taxa(core_list, phy_object)
    sample_core <- microbiome::transform(sample_core, transform = "compositional")

    sample_core %>%
        phyloseq::tax_glom(taxrank = "Genus") %>%
        phyloseq::psmelt() %>%
        tibble::as_tibble() %>%
        dplyr::mutate(across(Kingdom:Genus,
                             ~ gsub("^[a-z]__", "", .x))) %>%
        dplyr::select({{ tax_rank }}, Abundance) %>%
        dplyr::group_by({{ tax_rank }}) %>%
        dplyr::summarise(sum_abundance = sum(Abundance), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(sum_abundance))
}


#' @title taxa_colid
#' @description Tabla ancha de abundancias por taxon y muestra
#' @param df_taxa Tabla de taxones (output de taxa_table)
#' @param tax_rank Rango taxonomico (unquoted)
#' @param metadata Dataframe de metadatos para unir
#' @return Tibble en formato ancho unido con metadatos
taxa_colid <- function(df_taxa, tax_rank, metadata) {
    wide_tbl <- df_taxa %>%
        tidyr::drop_na() %>%
        dplyr::select(key, sample, Abundance, {{ tax_rank }}) %>%
        dplyr::group_by(key, {{ tax_rank }}) %>%
        dplyr::summarise(sum_abundance = sum(Abundance), .groups = "drop") %>%
        tidyr::pivot_wider(
            names_from = {{ tax_rank }},
            values_from = sum_abundance
        )
    dplyr::inner_join(metadata, wide_tbl, by = "key")
}

#' @title adiv_table
#' @description Tabla con métricas de diversidad alfa por muestra
#' @param phy_object Objeto phyloseq rarefaccionado con árbol
#' @param metric Vector de métricas para estimate_richness
#' @return Tibble con métricas de riqueza, diversidad, equitatividad y Faith PD
adiv_table <- function(phy_object, metric) {

    sample_data(phy_object)$no_reads <- sample_sums(phy_object)

    ## Riqueza y diversidad
    phy_adiv <- phy_object %>%
        phyloseq::estimate_richness(measures = metric) %>%
        tibble::rownames_to_column(var = "key")

    ## Faith PD con picante
    otu_mat <- as.data.frame(t(otu_table(phy_object)))
    tree    <- phy_tree(phy_object)
    faith   <- picante::pd(otu_mat, tree, include.root = FALSE) %>%
        tibble::rownames_to_column(var = "key") %>%
        dplyr::select(key, PD)

    ## Equitatividad
    phy_evenness <- phy_object %>%
        microbiome::evenness(index = "all") %>%
        tibble::rownames_to_column(var = "key")

    ## Metadatos
    data.frame(sample_data(phy_object)) %>%
        tibble::rownames_to_column(var = "key") %>%
        dplyr::full_join(phy_adiv,     by = "key") %>%
        dplyr::full_join(faith,        by = "key") %>%
        dplyr::full_join(phy_evenness, by = "key") %>%
        dplyr::mutate(dplyr::across(where(is.numeric), round, 2))
}
