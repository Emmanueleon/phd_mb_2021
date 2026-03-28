#' @title: Analisis diferencial con LEfSe
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Analisis diferencial de taxones usando LEfSe para identificar
#'   biomarcadores asociados a variables clinicas. Corre las combinaciones de
#'   los tres tipos de muestra y cinco variables clinicas, exportando un modelo
#'   y una tabla de biomarcadores por combinacion.

# Entradas:
#   out/phyloseq/phy_filtered.RDS

# Salidas (en out/tables/lefse/):
#   lefse_{sample}_{variable}_model.RDS
#   lefse_{sample}_{variable}_table.RDS
#   lefse_{sample}_{variable}_table.csv
#   lefse_summary.RDS
#   lefse_summary.csv

# Anterior: 09_diversidad_beta.R
# Siguiente: Pendiente

# ::::::::::::::::::::::::::::
# Preparacion del ambiente
# :::::::::::::::::::::::::::
rm(list = ls())

suppressPackageStartupMessages({
    library(phyloseq)
    library(tidyverse)
    library(microbiomeMarker)
    library(here)
})

## Crear directorio de salida
dir.create(here("out", "tables", "lefse"), recursive = TRUE, showWarnings = FALSE)

## Ruta de salida
out_dir <- here("out", "tables", "lefse")

## Cargar datos
message("Cargando datos...")
phy_filtered <- readRDS(here("out", "phyloseq", "phy_filtered.RDS"))
message("phy_filtered: ", nsamples(phy_filtered), " muestras")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Parametros del analisis
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
sample_types <- c("mom", "milk", "newborn")
variables <- c("imc", "fat", "delivery", "sex", "mom_age30")
taxa_rank <- "Genus"
sample_min_fixed <- 10

combinaciones <- tidyr::expand_grid(
    sample_type = sample_types,
    variable = variables
)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Funciones auxiliares
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
clean_lefse_table <- function(lefse_model) {
    marker_table(lefse_model) %>%
        as_tibble() %>%
        mutate(
            feature = str_replace(feature, "^[A-Za-z]_{1,2}", ""),
            feature = str_replace_all(feature, "_", " "),
            feature = str_squish(feature)
        ) %>%
        mutate(across(c(ef_lda, pvalue, padj), as.numeric)) %>%
        mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
        arrange(enrich_group, desc(ef_lda), feature) %>%
        mutate(order = row_number())
}

run_lefse_combo <- function(phy, sample_type_val, variable, taxa_rank = "Genus") {
    message("\n--- LEfSe: ", sample_type_val, " ~ ", variable, " ---")

    meta_df <- as.data.frame(sample_data(phy))
    if (!variable %in% colnames(meta_df)) {
        stop(
            "La variable '", variable, "' no existe en sample_data(phy)."
        )
    }

    keep_sample <- meta_df$sample == sample_type_val
    keep_value <- !is.na(meta_df[[variable]]) & meta_df[[variable]] != "na"
    phy_sub <- prune_samples(keep_sample & keep_value, phy)

    meta_sub <- as.data.frame(sample_data(phy_sub))
    group_counts <- table(meta_sub[[variable]])
    min_group_n <- if (length(group_counts) > 0) min(group_counts) else 0
    meets_sample_min <- min_group_n >= sample_min_fixed

    message("Muestras evaluadas: ", nsamples(phy_sub))
    message("Distribucion de grupos:")
    print(group_counts)

    if (nsamples(phy_sub) == 0) {
        stop("No hay muestras disponibles tras filtrar.")
    }

    if (length(group_counts) < 2) {
        stop("La variable no tiene al menos dos grupos tras filtrar.")
    }

    if (any(group_counts < 3)) {
        warning(
            "Hay grupos con n < 3 en ", sample_type_val, " ~ ", variable,
            ". Interpretar con cautela."
        )
    }

    if (!meets_sample_min) {
        stop(
            "La combinacion no cumple sample_min = ", sample_min_fixed,
            " (min_group_n = ", min_group_n, ")."
        )
    }

    lefse_model <- run_lefse(
        ps = phy_sub,
        group = variable,
        taxa_rank = taxa_rank,
        norm = "CPM",
        sample_min = sample_min_fixed,
        wilcoxon_cutoff = 0.05,
        lda_cutoff = 2
    )

    lefse_table <- clean_lefse_table(lefse_model)

    message("Biomarcadores detectados: ", nrow(lefse_table))

    list(
        sample_type = sample_type_val,
        variable = variable,
        taxa_rank = taxa_rank,
        n = nsamples(phy_sub),
        n_groups = length(group_counts),
        min_group_n = min_group_n,
        sample_min = sample_min_fixed,
        meets_sample_min = meets_sample_min,
        group_counts = group_counts,
        model = lefse_model,
        table = lefse_table
    )
}


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Analisis diferencial
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Analisis LEfSe para todas las combinaciones ---")

results <- pmap(
    combinaciones,
    function(sample_type, variable) {
        tryCatch(
            run_lefse_combo(
                phy = phy_filtered,
                sample_type_val = sample_type,
                variable = variable,
                taxa_rank = taxa_rank
            ),
            error = function(e) {
                message("  ERROR en ", sample_type, " ~ ", variable, ": ", e$message)
                list(
                    sample_type = sample_type,
                    variable = variable,
                    taxa_rank = taxa_rank,
                    n = NA_integer_,
                    n_groups = NA_integer_,
                    min_group_n = NA_integer_,
                    sample_min = sample_min_fixed,
                    meets_sample_min = NA,
                    group_counts = NA,
                    model = NULL,
                    table = tibble(),
                    error = e$message
                )
            }
        )
    }
)

result_names <- paste(combinaciones$sample_type, combinaciones$variable, sep = "_")
results <- set_names(results, result_names)

lefse_summary <- imap_dfr(results, function(r, result_name) {
    group_counts_chr <- if (is.table(r$group_counts)) {
        paste(names(r$group_counts), as.integer(r$group_counts), collapse = ", ")
    } else {
        NA_character_
    }

    tibble(
        analysis_id = result_name,
        sample_type = r$sample_type,
        variable = r$variable,
        taxa_rank = r$taxa_rank,
        n = r$n,
        n_groups = r$n_groups,
        min_group_n = r$min_group_n,
        sample_min = r$sample_min,
        meets_sample_min = r$meets_sample_min,
        biomarkers = nrow(r$table),
        group_counts = group_counts_chr,
        status = ifelse(is.null(r$model), "error", "ok"),
        error = ifelse(is.null(r$model), r$error, NA_character_)
    )
})

print(lefse_summary)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Exportacion
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Exportacion ---")

walk(results, function(r) {
    if (is.null(r$model)) {
        return()
    }

    base_name <- paste("lefse", r$sample_type, r$variable, sep = "_")

    saveRDS(
        r$model,
        file.path(out_dir, paste0(base_name, "_model.RDS"))
    )
    saveRDS(
        r$table,
        file.path(out_dir, paste0(base_name, "_table.RDS"))
    )
    write.csv(
        r$table,
        file.path(out_dir, paste0(base_name, "_table.csv")),
        row.names = FALSE
    )

    message("  ", base_name, " exportado")
})

saveRDS(lefse_summary, file.path(out_dir, "lefse_summary.RDS"))
write.csv(lefse_summary, file.path(out_dir, "lefse_summary.csv"), row.names = FALSE)
message("  lefse_summary exportado")

message("Archivos en: ", out_dir)
message("\n--- 10_differential_analysis.R completado ---")
