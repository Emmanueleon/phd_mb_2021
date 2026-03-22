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
    otu <- phyloseq::otu_table(as(biom_data, "matrix"), taxa_are_rows = TRUE)

    ## Taxonomía — parsear niveles separados por ";"
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
            sep = ";\\s*",
            fill = "right"
        ) %>%
        tibble::column_to_rownames("Feature ID") %>%
        dplyr::select(-Confidence) %>%
        as.matrix()

    tax <- phyloseq::tax_table(tax_mat)

    ## Árbol filogenético
    tree <- ape::read.tree(file.path(export_dir, "tree.nwk"))

    ## Metadatos
    meta_raw <- readr::read_tsv(
        metadata_path,
        show_col_types = FALSE,
        comment        = "#"
    ) %>%
        dplyr::rename(sample_id = `sample-id`) %>%
        tibble::column_to_rownames("sample_id")

    meta <- phyloseq::sample_data(meta_raw)

    ## Construir phyloseq
    ps <- phyloseq::phyloseq(otu, tax, tree, meta)
    return(ps)
}


# =============================================================================
# 2. ESTADÍSTICA GENERAL
# =============================================================================

#' desc_stats
#' @description Calcula estadística descriptiva de un dataframe numérico
#' @param df Dataframe con variables numéricas (usar select() previamente)
#' @return Tibble con n, media, sd, mediana, min, max, Q1, SE, IC95
desc_stats <- function(df) {
    df %>%
        purrr::map_df(~ tibble::tibble(
            no        = length(.x),
            mean      = mean(.x, na.rm = TRUE),
            sd        = sd(.x, na.rm = TRUE),
            median    = median(.x, na.rm = TRUE),
            min       = min(.x, na.rm = TRUE),
            max       = max(.x, na.rm = TRUE),
            q25       = quantile(.x, probs = 0.25, na.rm = TRUE),
            q75       = quantile(.x, probs = 0.75, na.rm = TRUE),
            se        = sd / sqrt(no),
            lower.ci  = mean - qt(1 - (0.05 / 2), no - 1) * se,
            upper.ci  = mean + qt(1 - (0.05 / 2), no - 1) * se
        ), .id = "variable") %>%
        dplyr::mutate(dplyr::across(where(is.numeric), round, 2))
}


#' descshap_tbl
#' @description Estadística descriptiva + prueba de Shapiro-Wilk por grupo
#' @param df Dataframe
#' @param var Variable de agrupación (unquoted)
#' @param no_group Número de grupo a evaluar (índice)
#' @param var_vector Vector de nombres de variables a evaluar
#' @return Tibble con descriptivos, p-valor Shapiro y distribución
descshap_tbl <- function(df, var, no_group, var_vector) {
    df_nested <- df %>%
        dplyr::group_by({{ var }}) %>%
        tidyr::nest() %>%
        dplyr::arrange({{ var }})

    df_plucked <- df_nested %>%
        purrr::pluck("data", no_group) %>%
        dplyr::select(dplyr::all_of(var_vector))

    shap_test <- df_plucked %>%
        purrr::map(shapiro.test) %>%
        purrr::map_df(broom::tidy, .id = "variable") %>%
        dplyr::mutate(distribution = ifelse(p.value >= 0.05, "Parametrica", "No Parametrica"))

    desc_tbl <- df_plucked %>%
        desc_stats() %>%
        dplyr::left_join(shap_test, by = "variable") %>%
        dplyr::mutate(group = df_nested[[no_group, 1]])

    return(desc_tbl)
}


#' inferential_tests
#' @description Comparación entre 2 grupos con corrección FDR
#' @param df Dataframe
#' @param var Vector de nombres de variables a comparar
#' @param method Función del test (ej. t.test, wilcox.test)
#' @return Tibble con estadístico, p-valor, p-ajustado y significancia
inferential_tests <- function(df, var, method) {
    df %>%
        dplyr::select(dplyr::all_of(var)) %>%
        purrr::map_df(~ broom::tidy(method(. ~ group_var)), .id = "variable") %>%
        rstatix::add_significance() %>%
        rstatix::adjust_pvalue(p.col = "p.value", method = "fdr") %>%
        rstatix::add_significance() %>%
        dplyr::mutate(dplyr::across(where(is.numeric), round, 5)) %>%
        dplyr::mutate(stars = dplyr::case_when(
            p.value < 0.001 ~ "***",
            dplyr::between(p.value, 0.001, 0.009) ~ "**",
            dplyr::between(p.value, 0.01, 0.049) ~ "*",
            dplyr::between(p.value, 0.05, 0.079) ~ "#",
            p.value > 0.08 ~ "ns",
            TRUE ~ ""
        ))
}


# =============================================================================
# 3. CORRELACIONES
# =============================================================================

#' cor_df
#' @description Tabla de correlación con r, n y p (Pearson o Spearman)
#' @param df Dataframe numérico
#' @param method "pearson" o "spearman"
#' @return Tibble en formato largo con r, P, n y flags de significancia
cor_df <- function(df, method = "pearson") {
    M <- Hmisc::rcorr(as.matrix(df), type = method)

    purrr::map(M, ~ data.frame(.x)) %>%
        purrr::map(~ tibble::rownames_to_column(.x, var = "measure1")) %>%
        purrr::map(~ tidyr::pivot_longer(.x, -measure1,
            names_to = "measure2",
            values_to = "value"
        )) %>%
        dplyr::bind_rows(.id = "id") %>%
        tidyr::pivot_wider(names_from = id, values_from = value) %>%
        dplyr::mutate(
            sig_p    = ifelse(P < 0.05, TRUE, FALSE),
            p_if_sig = ifelse(P < 0.05, P, NA),
            r_if_sig = ifelse(P < 0.05, r, NA)
        )
}

# Wrappers para compatibilidad con código anterior
#' @describeIn cor_df Correlación de Pearson (wrapper)
cors_pearson <- function(df) cor_df(df, method = "pearson")

#' @describeIn cor_df Correlación de Spearman (wrapper)
cors_spearman <- function(df) cor_df(df, method = "spearman")

#' @describeIn cor_df Tabla formateada de correlación (Pearson o Spearman)
formatted_cor <- function(df, method = "pearson") cor_df(df, method = method)

# Alias para compatibilidad con código anterior
formatted_corpea <- function(df) cor_df(df, method = "pearson")
formatted_corsper <- function(df) cor_df(df, method = "spearman")


# =============================================================================
# 4. UTILIDADES GRÁFICAS
# =============================================================================

#' get_only_legend
#' @description Extrae la leyenda de un gráfico ggplot
#' @param plot Objeto ggplot
#' @return Grob con la leyenda
get_only_legend <- function(plot) {
    plot_table <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot))
    legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
    plot_table$grobs[[legend_plot]]
}


# =============================================================================
# 5. CALIDAD Y FILTRADO
# =============================================================================

#' reads_tbl
#' @description Tabla con número de lecturas por muestra
#' @param phy_object Objeto phyloseq
#' @param sample Vector con los niveles de sample_type
#' @return Tibble con sample_id, sample_type y reads
reads_tbl <- function(phy_object, sample) {
    phyloseq::sample_data(phy_object)$reads <- phyloseq::sample_sums(phy_object)
    as.data.frame(phyloseq::sample_data(phy_object)) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(sample_type = factor(sample_type, levels = sample)) %>%
        dplyr::select(sample_id, sample_type, reads)
}


#' histogram_plot
#' @description Histograma de distribución de lecturas por tipo de muestra
#' @param phy_object Objeto phyloseq
#' @param breaks Vector de niveles de sample_type
#' @param color Vector de colores
#' @param labels Etiquetas de los grupos
#' @param title Título del gráfico
#' @param subtitle Subtítulo del gráfico
#' @param labeller Named vector para facet labels
#' @return Objeto ggplot
histogram_plot <- function(phy_object, breaks, color, labels, title, subtitle, labeller) {
    phyloseq::sample_data(phy_object)$reads <- phyloseq::sample_sums(phy_object)

    as.data.frame(phyloseq::sample_data(phy_object)) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(sample_type = factor(sample_type, levels = breaks)) %>%
        ggplot2::ggplot(ggplot2::aes(x = reads, fill = sample_type)) +
        ggplot2::geom_histogram(
            alpha = 0.8, binwidth = 10000,
            breaks = seq(0, 90000, by = 5000)
        ) +
        ggplot2::geom_vline(
            xintercept = 5000, linetype = "dashed",
            color = "red", linewidth = 1.2
        ) +
        ggplot2::scale_y_continuous(
            expand = c(0, 0), limits = c(0, 10),
            breaks = seq(0, 10, by = 2)
        ) +
        ggplot2::scale_x_continuous(
            expand = c(0, 0), limits = c(0, 90000),
            labels = scales::label_number(),
            breaks = seq(0, 90000, by = 10000)
        ) +
        ggplot2::scale_fill_manual(
            name = NULL, breaks = breaks,
            values = color, labels = labels
        ) +
        ggplot2::labs(
            title = title, subtitle = subtitle,
            x = "\nTamaño de la librería (lecturas)",
            y = "No. muestras\n"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::facet_wrap(~sample_type,
            ncol = 1,
            labeller = ggplot2::as_labeller(labeller)
        ) +
        ggplot2::theme(
            text = ggplot2::element_text(size = 22, color = "black"),
            legend.position = "none",
            plot.title = ggplot2::element_text(size = 20, face = "bold"),
            strip.background = ggplot2::element_rect(color = "gray", fill = "gray"),
            strip.text.x = ggplot2::element_text(size = 15, color = "black"),
            plot.title.position = "plot"
        )
}


#' rarecurve_plot
#' @description Curva de rarefacción para evaluar profundidad de secuenciación
#' @param phy_object Objeto phyloseq
#' @param var Variable de agrupación (string)
#' @param breaks Vector de niveles
#' @param var_color Vector de colores
#' @param labels Etiquetas
#' @param title Título
#' @param subtitle Subtítulo
#' @return Objeto ggplot
rarecurve_plot <- function(phy_object, var, breaks, var_color, labels, title, subtitle) {
    graph_rarcurve <- phyloseq.extended::ggrare(phy_object,
        step = 50,
        color = var, label = NULL, se = FALSE
    )
    graph_rarcurve +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 600)) +
        ggplot2::scale_x_continuous(
            expand = c(0, 0), limits = c(0, 90000),
            labels = scales::label_number(),
            breaks = seq(0, 90000, by = 10000)
        ) +
        ggplot2::scale_colour_manual(
            name = NULL, breaks = breaks,
            values = var_color, labels = labels
        ) +
        ggplot2::labs(
            title = title, subtitle = subtitle,
            x = "\nTamaño de la librería (lecturas)",
            y = "Riqueza (No. de especies)\n"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            text = ggplot2::element_text(size = 15, color = "black"),
            plot.title = ggplot2::element_text(size = 17, face = "bold"),
            plot.title.position = "plot"
        )
}


# =============================================================================
# 6. DIVERSIDAD ALFA
# =============================================================================

#' adiv_table
#' @description Tabla con métricas de diversidad alfa por muestra
#' @param phy_object Objeto phyloseq (rarefaccionado)
#' @param metric Vector de métricas (ej. c("Chao1", "Shannon"))
#' @return Tibble con métricas de riqueza, diversidad, equitatividad y Faith PD
adiv_table <- function(phy_object, metric) {
    phyloseq::sample_data(phy_object)$no_reads <- phyloseq::sample_sums(phy_object)

    phy_adiv <- phy_object %>%
        phyloseq::estimate_richness(measures = metric) %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::mutate(microbiome::rarity(phy_object))

    faith_pd <- phy_object %>%
        btools::estimate_pd() %>%
        tibble::rownames_to_column(var = "sample_id")

    phy_evenness <- phy_object %>%
        microbiome::evenness(index = "all") %>%
        tibble::rownames_to_column(var = "sample_id")

    as.data.frame(phyloseq::sample_data(phy_object)) %>%
        tibble::as_tibble() %>%
        dplyr::select(sample_id, no_reads) %>%
        dplyr::full_join(phy_adiv, by = "sample_id") %>%
        dplyr::full_join(faith_pd, by = "sample_id") %>%
        dplyr::full_join(phy_evenness, by = "sample_id") %>%
        dplyr::mutate(dplyr::across(where(is.numeric), round, 2)) %>%
        dplyr::mutate(sample_id = factor(sample_id))
}


#' adiv_distribution
#' @description Distribución de métricas alfa por tipo de muestra y variable
#' @param df_adiv Tabla de diversidad alfa (output de adiv_table)
#' @param var Variable de agrupación (unquoted)
#' @param metric Vector de métricas a evaluar
#' @return Tibble con distribución y prueba de Shapiro-Wilk
adiv_distribution <- function(df_adiv, var, metric) {
    df_adiv %>%
        dplyr::select(sample_type, {{ var }}, dplyr::all_of(metric)) %>%
        tidyr::pivot_longer(dplyr::all_of(metric),
            names_to = "metric", values_to = "value"
        ) %>%
        dplyr::filter(!is.na(value), !is.na({{ var }})) %>%
        dplyr::group_by(sample_type, {{ var }}, metric) %>%
        tidyr::nest() %>%
        dplyr::mutate(distribution = purrr::map(
            data, ~ shapiro.test(.x$value) %>% broom::tidy()
        )) %>%
        tidyr::unnest(distribution) %>%
        dplyr::mutate(distribucion = ifelse(p.value > 0.05, "Parametrica", "No Parametrica")) %>%
        dplyr::arrange(sample_type, metric, {{ var }}) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(dplyr::across(where(is.numeric), round, 3))
}


#' ametric_plot
#' @description Gráfica de métrica de diversidad alfa por tipo de muestra y variable
#' @param df Tabla larga de métricas alfa
#' @param metric2graph Vector de métricas a graficar
#' @param var Variable de agrupación (unquoted)
#' @param method Método de comparación (ej. "wilcox.test")
#' @param color Vector de colores
#' @param breaks Vector de niveles
#' @param labels Etiquetas
#' @param title Título
#' @param subtitle Subtítulo
#' @param index Nombre del eje Y
#' @param labeller Named vector para facet labels
#' @return Objeto ggplot
ametric_plot <- function(df, metric2graph, var, method, color, breaks,
                         labels, title, subtitle, index, labeller) {
    df %>%
        dplyr::filter(metric %in% metric2graph) %>%
        ggplot2::ggplot(ggplot2::aes(x = {{ var }}, y = value)) +
        ggplot2::stat_summary(
            fun.data = ggplot2::median_hilow, fun.args = 0.50,
            geom = "crossbar", width = 0.5, alpha = 0.5,
            show.legend = FALSE
        ) +
        ggpubr::stat_compare_means(
            method = method, size = 3,
            label.x.npc = .9, label.y.npc = 0.99
        ) +
        ggplot2::geom_jitter(ggplot2::aes(fill = {{ var }}),
            height = 0, width = 0.25, shape = 21,
            size = 3, alpha = 0.8, color = "black",
            show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(
            name = NULL, values = color,
            breaks = breaks, labels = labels
        ) +
        ggplot2::scale_x_discrete(breaks = breaks, labels = labels) +
        ggplot2::labs(title = title, subtitle = subtitle, x = "", y = index) +
        ggplot2::facet_wrap(~sample_type,
            labeller = ggplot2::as_labeller(labeller)
        ) +
        cowplot::theme_cowplot() +
        ggplot2::theme(
            text = ggplot2::element_text(size = 15, color = "black"),
            legend.position = "none",
            plot.title = ggplot2::element_text(size = 17, face = "bold"),
            strip.background = ggplot2::element_rect(color = "gray", fill = "gray"),
            strip.text.x = ggplot2::element_text(size = 15, color = "black"),
            plot.title.position = "plot"
        )
}


#' asample_plot
#' @description Gráfica de métrica alfa para una sola muestra con comparaciones
#' @param df Tabla de métricas alfa filtrada por muestra
#' @param metric2graph Métrica a graficar (unquoted)
#' @param var Variable de agrupación (unquoted)
#' @param method Método de comparación
#' @param my_comparisons Lista de comparaciones
#' @param ylimit Límite superior del eje Y
#' @param ylabel Posición Y de las etiquetas de comparación
#' @param color Vector de colores
#' @param breaks Vector de niveles
#' @param labels Etiquetas
#' @param title Título
#' @param subtitle Subtítulo
#' @param index Nombre del eje Y
#' @return Objeto ggplot
asample_plot <- function(df, metric2graph, var, method, my_comparisons,
                         ylimit, ylabel, color, breaks, labels,
                         title, subtitle, index) {
    df %>%
        ggplot2::ggplot(ggplot2::aes(
            x = {{ var }}, y = {{ metric2graph }},
            fill = {{ var }}
        )) +
        ggplot2::stat_summary(
            fun.data = ggplot2::median_hilow, fun.args = 0.50,
            geom = "crossbar", width = 0.5, alpha = 0.5,
            show.legend = FALSE
        ) +
        ggpubr::stat_compare_means(
            method = method, comparisons = my_comparisons,
            label = "p.signif", bracket.size = 0.7,
            color = "black", size = 5, label.y = ylabel
        ) +
        ggplot2::geom_jitter(ggplot2::aes(fill = {{ var }}),
            height = 0, width = 0.25, shape = 21,
            size = 3, alpha = 0.8, color = "black",
            show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(
            name = NULL, values = color,
            breaks = breaks, labels = labels
        ) +
        ggplot2::scale_x_discrete(breaks = breaks, labels = labels) +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, ylimit)) +
        ggplot2::labs(title = title, subtitle = subtitle, x = "", y = index) +
        cowplot::theme_cowplot() +
        ggplot2::theme(
            text = ggplot2::element_text(size = 15, color = "black"),
            legend.position = "none",
            plot.title = ggplot2::element_text(size = 17, face = "bold"),
            plot.title.position = "plot"
        )
}


#' ametrics_plots
#' @description Gráfica de múltiples métricas alfa para una sola muestra
#' @param df Tabla larga de métricas alfa
#' @param sample Tipo de muestra a evaluar
#' @param var Variable de agrupación (unquoted)
#' @param method Método de comparación
#' @param color Vector de colores
#' @param breaks Vector de niveles
#' @param labels Etiquetas
#' @param title Título
#' @param subtitle Subtítulo
#' @return Objeto ggplot
ametrics_plots <- function(df, sample, var, method, color, breaks,
                           labels, title, subtitle) {
    df %>%
        dplyr::filter(sample_type == sample) %>%
        ggplot2::ggplot(ggplot2::aes(x = {{ var }}, y = value)) +
        ggplot2::stat_summary(
            fun.data = ggplot2::median_hilow, fun.args = 0.50,
            geom = "crossbar", width = 0.5, alpha = 0.5,
            show.legend = FALSE
        ) +
        ggpubr::stat_compare_means(
            method = method, size = 3,
            label.x.npc = .9, label.y.npc = 0.99
        ) +
        ggplot2::geom_jitter(ggplot2::aes(fill = {{ var }}),
            height = 0, width = 0.25, shape = 21,
            size = 3, alpha = 0.8, color = "black",
            show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(
            name = NULL, values = color,
            breaks = breaks, labels = labels
        ) +
        ggplot2::scale_x_discrete(breaks = breaks, labels = labels) +
        ggplot2::labs(title = title, subtitle = subtitle, x = "", y = "") +
        ggplot2::facet_wrap(~metric, scales = "free") +
        ggplot2::theme(
            text = ggplot2::element_text(size = 15, color = "black"),
            legend.position = "none",
            plot.title = ggplot2::element_text(size = 17, face = "bold"),
            strip.background = ggplot2::element_rect(color = "gray"),
            strip.text.x = ggplot2::element_text(size = 15, color = "black"),
            plot.title.position = "plot"
        )
}


# =============================================================================
# 7. DIVERSIDAD BETA
# =============================================================================

#' bdiv_plot
#' @description Gráfica de ordenación para diversidad beta (PCoA, NMDS, etc.)
#' @param phy_object Objeto phyloseq
#' @param ordination Método de ordenación ("PCoA", "NMDS", etc.)
#' @param bdistance Distancia ("unifrac", "wunifrac", "bray", etc.)
#' @param var Variable de agrupación (string)
#' @param breaks Vector de niveles
#' @param values Vector de colores
#' @param labels Etiquetas
#' @return Objeto ggplot
bdiv_plot <- function(phy_object, ordination, bdistance, var,
                      breaks, values, labels) {
    ord <- phyloseq::ordinate(phy_object,
        method = ordination,
        distance = bdistance
    )
    p <- phyloseq::plot_ordination(phy_object, ord, col = var)

    p +
        ggplot2::geom_point(size = 2, alpha = 0.75) +
        ggplot2::stat_ellipse(show.legend = FALSE) +
        ggplot2::scale_color_manual(
            name = NULL, breaks = breaks,
            values = values, labels = labels
        ) +
        cowplot::theme_cowplot() +
        ggplot2::theme(
            text = ggplot2::element_text(size = 22, color = "black"),
            legend.position = "right",
            plot.title = ggplot2::element_text(size = 22, face = "bold"),
            plot.title.position = "plot"
        )
}


#' p.adjust.envfit
#' @description Ajusta los p-valores de un objeto envfit
#' @author David Zelený
#' @param x Objeto envfit
#' @param method Método de ajuste (default: "bonferroni")
#' @param n Número total de tests (opcional)
#' @return Objeto envfit con p-valores ajustados
p.adjust.envfit <- function(x, method = "bonferroni", n) {
    x.new <- x
    pval.vectors <- if (!is.null(x$vectors)) x$vectors$pvals else NULL
    pval.factors <- if (!is.null(x$factors)) x$factors$pvals else NULL
    if (missing(n)) n <- length(pval.vectors) + length(pval.factors)
    if (!is.null(x$vectors)) x.new$vectors$pvals <- p.adjust(x$vectors$pvals, method, n)
    if (!is.null(x$factors)) x.new$factors$pvals <- p.adjust(x$factors$pvals, method, n)
    cat("Adjustment of significance by", method, "method\n")
    return(x.new)
}


# =============================================================================
# 8. COMPOSICIÓN TAXONÓMICA
# =============================================================================

#' taxa_table
#' @description Crea tabla larga con abundancias y rangos taxonómicos por muestra
#' @param phy_object Objeto phyloseq
#' @return Tibble con sample_id, Abundance, OTU, Kingdom:Species
taxa_table <- function(phy_object) {
    phy_object %>%
        phyloseq::psmelt() %>%
        tibble::as_tibble() %>%
        dplyr::select(sample_id, Abundance, OTU, Kingdom:Species) %>%
        dplyr::mutate(Kingdom = gsub("d__", "", Kingdom)) %>%
        tidyr::replace_na(list(
            Kingdom = "Unclassified", Phylum  = "Unclassified",
            Class   = "Unclassified", Order   = "Unclassified",
            Family  = "Unclassified", Genus   = "Unclassified",
            Species = "Unclassified"
        ))
}


#' top_taxa_lst
#' @description Identifica los top n taxones más abundantes
#' @param df_taxa Tabla de taxones (output de taxa_table)
#' @param tax_rank Rango taxonómico (unquoted, ej. Genus)
#' @param top Número de taxones a identificar
#' @return Vector con los nombres de los top taxones
top_taxa_lst <- function(df_taxa, tax_rank, top) {
    df_taxa %>%
        dplyr::select({{ tax_rank }}, Abundance) %>%
        dplyr::group_by({{ tax_rank }}) %>%
        dplyr::summarise(sum_abun = sum(Abundance)) %>%
        dplyr::arrange(dplyr::desc(sum_abun)) %>%
        dplyr::distinct({{ tax_rank }}) %>%
        dplyr::slice_head(n = top) %>%
        dplyr::pull()
}


#' taxa_lst
#' @description Identifica los top n taxones por filo
#' @param df_taxa Tabla de taxones
#' @param tax_rank Rango taxonómico (unquoted)
#' @param top Número de taxones por filo
#' @param phylum_vector Vector de filos de interés
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


#' taxa_colid
#' @description Tabla ancha de abundancias por taxón y muestra
#' @param df_taxa Tabla de taxones (output de taxa_table)
#' @param tax_rank Rango taxonómico (unquoted)
#' @param metadata Dataframe de metadatos para unir (evita dependencia global)
#' @return Tibble en formato ancho unido con metadatos
taxa_colid <- function(df_taxa, tax_rank, metadata) {
    wide_tbl <- df_taxa %>%
        tidyr::drop_na() %>%
        dplyr::select(sample_id, Abundance, {{ tax_rank }}) %>%
        dplyr::group_by(sample_id, {{ tax_rank }}) %>%
        dplyr::summarise(sum_abundance = sum(Abundance), .groups = "drop") %>%
        tidyr::pivot_wider(
            names_from = {{ tax_rank }},
            values_from = sum_abundance
        )
    dplyr::inner_join(metadata, wide_tbl, by = "sample_id")
}


#' taxa_summ
#' @description Igual que taxa_colid pero con metadatos unidos
#' @note Wrapper de taxa_colid para compatibilidad con código anterior
#' @param df_taxa Tabla de taxones
#' @param tax_rank Rango taxonómico (unquoted)
#' @param metadata Dataframe de metadatos
#' @return Tibble en formato ancho con metadatos
taxa_summ <- function(df_taxa, tax_rank, metadata) {
    taxa_colid(df_taxa, {{ tax_rank }}, metadata)
}


#' taxbar_by_id
#' @description Gráfica de barras de composición por muestra individual
#' @param df Tabla larga con taxones y abundancias
#' @param taxa_fill Rango taxonómico para colorear (unquoted)
#' @param taxa_colors Vector de colores
#' @return Objeto ggplot
taxbar_by_id <- function(df, taxa_fill, taxa_colors) {
    order_taxa_id <- df %>%
        dplyr::filter(prevalence_phylum == "Firmicutes") %>%
        dplyr::count(sample_id, wt = Abundance, sort = TRUE) %>%
        dplyr::pull(sample_id) %>%
        as.character()

    df %>%
        dplyr::mutate(dplyr::across(
            Phylum:prevalence_genus,
            ~ stringr::str_replace_all(.x, "_", " ")
        )) %>%
        dplyr::mutate(sample_id = forcats::fct_relevel(sample_id, order_taxa_id)) %>%
        ggplot2::ggplot(ggplot2::aes(
            x = sample_id, y = Abundance * 100,
            fill = {{ taxa_fill }}
        )) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
        ggplot2::scale_fill_manual(name = NULL, values = taxa_colors)
}


#' taxbar_by_group
#' @description Gráfica de barras de composición por variable de agrupación
#' @param df Tabla larga con taxones
#' @param var Variable de agrupación (unquoted)
#' @param levels Niveles de la variable
#' @param labels Etiquetas
#' @param phylum_order Filo para ordenar muestras
#' @param taxa_fill Rango taxonómico para colorear (unquoted)
#' @param taxa_colors Vector de colores
#' @param title Título
#' @param subtitle Subtítulo
#' @return Objeto ggplot
taxbar_by_group <- function(df, var, levels, labels, phylum_order,
                            taxa_fill, taxa_colors, title, subtitle) {
    order_taxa <- df %>%
        dplyr::select(sample_id, Abundance, prevalence_phylum:prevalence_genus) %>%
        dplyr::group_by(sample_id, prevalence_phylum) %>%
        dplyr::summarise(mean_abundance = sum(Abundance), .groups = "drop") %>%
        dplyr::filter(prevalence_phylum == phylum_order) %>%
        dplyr::arrange(dplyr::desc(mean_abundance)) %>%
        dplyr::mutate(order = dplyr::row_number()) %>%
        dplyr::select(sample_id, order)

    df %>%
        dplyr::inner_join(order_taxa, by = "sample_id") %>%
        tidyr::drop_na(Abundance) %>%
        dplyr::mutate(
            sample_id = forcats::fct_reorder(factor(sample_id), order),
            var       = factor({{ var }}, levels = levels, labels = labels)
        ) %>%
        ggplot2::ggplot(ggplot2::aes(
            x = sample_id, y = Abundance * 100,
            fill = {{ taxa_fill }}
        )) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
        ggplot2::scale_fill_manual(name = NULL, values = taxa_colors) +
        ggplot2::labs(
            title = title, subtitle = subtitle,
            fill = "Phylum", y = "Abundancia Relativa (%)\n"
        ) +
        ggplot2::facet_wrap(~var, scales = "free_x", nrow = 2) +
        cowplot::theme_cowplot() +
        ggplot2::theme(
            legend.position = "right",
            plot.title = ggplot2::element_text(size = 17, face = "bold"),
            strip.background = ggplot2::element_rect(color = "gray"),
            strip.text.x = ggplot2::element_text(size = 12, color = "black"),
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 9, angle = 90),
            plot.title.position = "plot"
        )
}


#' taxbar_by_relabun
#' @description Gráfica de barras de abundancia relativa promedio por grupo
#' @param df Tabla larga con taxones
#' @param var Variable de agrupación (unquoted)
#' @param taxa_fill Rango taxonómico para colorear (unquoted)
#' @param breaks Niveles de la variable
#' @param levels Etiquetas de los niveles
#' @param taxa_colors Vector de colores
#' @param title Título
#' @param subtitle Subtítulo
#' @return Objeto ggplot
taxbar_by_relabun <- function(df, var, taxa_fill, breaks, levels,
                              taxa_colors, title, subtitle) {
    samtaxa_sum <- df %>%
        dplyr::group_by({{ var }}, sample_id, Phylum, Family, Genus, {{ taxa_fill }}) %>%
        dplyr::summarise(rel_abund = sum(Abundance), .groups = "drop") %>%
        dplyr::group_by({{ var }}, Phylum, Family, Genus, {{ taxa_fill }}) %>%
        dplyr::summarise(mean_rel_abun = 100 * mean(rel_abund), .groups = "drop")

    samtaxa_sum %>%
        tidyr::drop_na(mean_rel_abun) %>%
        dplyr::mutate(var = factor({{ var }}, levels = breaks, labels = levels)) %>%
        ggplot2::ggplot(ggplot2::aes(
            x = var, y = mean_rel_abun,
            fill = {{ taxa_fill }}
        )) +
        ggplot2::geom_col() +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::scale_fill_manual(name = NULL, values = taxa_colors) +
        ggplot2::labs(
            title = title, subtitle = subtitle,
            x = NULL, y = "Abundancia relativa (%)"
        ) +
        cowplot::theme_cowplot() +
        ggplot2::theme(
            legend.position = "right",
            plot.title = ggplot2::element_text(size = 17, face = "bold"),
            plot.title.position = "plot"
        )
}


#' taxa2_compare_plots
#' @description Boxplots de abundancia de taxones por grupo
#' @param df Tabla con taxones
#' @param sample Tipo de muestra
#' @param var Variable de agrupación (unquoted)
#' @param color Vector de colores
#' @param levels Niveles
#' @param breaks Breaks
#' @param labels Etiquetas
#' @param title Título
#' @param subtitle Subtítulo
#' @param taxa Variable de taxones para facet
#' @return Objeto ggplot
taxa2_compare_plots <- function(df, sample, var, color, levels, breaks,
                                labels, title, subtitle, taxa) {
    df %>%
        dplyr::filter(sample_type == sample) %>%
        dplyr::mutate(var = factor({{ var }}, levels = levels)) %>%
        ggplot2::ggplot(ggplot2::aes(x = var, y = Abundance * 100)) +
        ggplot2::stat_summary(
            fun.data = ggplot2::median_hilow, fun.args = 0.50,
            geom = "crossbar", width = 0.5, alpha = 0.5,
            show.legend = FALSE
        ) +
        ggpubr::stat_compare_means(
            method = "wilcox.test", size = 3,
            label.x.npc = .9, label.y.npc = 0.99
        ) +
        ggplot2::geom_jitter(ggplot2::aes(fill = {{ var }}),
            height = 0, width = 0.25, shape = 21,
            size = 3, alpha = 0.8, color = "black",
            show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(
            name = NULL, values = color,
            breaks = breaks, labels = labels
        ) +
        ggplot2::scale_x_discrete(breaks = breaks, labels = labels) +
        ggplot2::labs(title = title, subtitle = subtitle, x = "", y = "") +
        ggplot2::facet_wrap(~taxa, scales = "free") +
        ggplot2::theme(
            text = ggplot2::element_text(size = 15, color = "black"),
            legend.position = "right",
            plot.title = ggplot2::element_text(size = 17, face = "bold"),
            strip.background = ggplot2::element_rect(color = "gray"),
            strip.text.x = ggplot2::element_text(size = 15, color = "black"),
            plot.title.position = "plot"
        )
}


# =============================================================================
# 9. CORE MICROBIOTA
# =============================================================================

#' core_taxa
#' @description Tabla de abundancia de taxones del core microbiota
#' @param phy_object Objeto phyloseq
#' @param preval Umbral de prevalencia (ej. 0.70 para 70%)
#' @param tax_rank Rango taxonómico (unquoted)
#' @return Tibble con taxones del core y su abundancia total
core_taxa <- function(phy_object, preval, tax_rank) {
    core_list <- microbiome::core_members(phy_object, prevalence = preval)
    sample_core <- phyloseq::prune_taxa(core_list, phy_object)
    sample_core <- microbiome::transform(sample_core, transform = "compositional")

    sample_core %>%
        phyloseq::tax_glom(taxrank = "Genus") %>%
        phyloseq::psmelt() %>%
        tibble::as_tibble() %>%
        dplyr::select({{ tax_rank }}, Abundance) %>%
        dplyr::group_by({{ tax_rank }}) %>%
        dplyr::summarise(sum_abundance = sum(Abundance), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(sum_abundance))
}


# =============================================================================
# 10. ABUNDANCIA DIFERENCIAL
# =============================================================================

#' taxa_relabun_by_group
#' @description Tabla de abundancia relativa promedio por grupo y rango taxonómico
#' @param df Tabla larga con taxones
#' @param var Variable de agrupación (unquoted)
#' @param taxa_fill Rango taxonómico para colorear (unquoted)
#' @param group_var Valor del grupo a filtrar (string)
#' @param taxa_rank Rango taxonómico para agrupar (unquoted)
#' @return Dos tablas impresas: abundancias promedio y suma por grupo
taxa_relabun_by_group <- function(df, var, taxa_fill, group_var, taxa_rank) {
    samtaxa_sum <- df %>%
        dplyr::group_by({{ var }}, sample_id, Phylum, Family, Genus, {{ taxa_fill }}) %>%
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


#' relabun_tbl
#' @description Tabla de abundancia relativa por rango taxonómico y grupo
#' @param df Tabla larga con taxones
#' @param var Variable de agrupación (unquoted)
#' @param var_group Valor del grupo a filtrar (string)
#' @param tax_rank Rango taxonómico (unquoted)
#' @return Tibble con abundancia relativa por taxón
relabun_tbl <- function(df, var, var_group, tax_rank) {
    taxon_tbl <- df %>%
        dplyr::group_by(
            {{ var }}, sample_id, Phylum, Family, Genus,
            prevalence_phylum, prevalence_family, prevalence_genus
        ) %>%
        dplyr::summarise(rel_abund = sum(Abundance), .groups = "drop") %>%
        dplyr::group_by(
            {{ var }}, Phylum, Family, Genus,
            prevalence_phylum, prevalence_family, prevalence_genus
        ) %>%
        dplyr::summarise(mean_rel_abun = 100 * mean(rel_abund), .groups = "drop") %>%
        dplyr::filter({{ var }} == var_group)

    taxon_tbl %>%
        dplyr::group_by({{ tax_rank }}) %>%
        dplyr::summarise(sum_abun = sum(mean_rel_abun)) %>%
        dplyr::arrange(dplyr::desc(sum_abun))
}


message(
    "functions_microbiota.R cargado correctamente — ",
    length(lsf.str()), " funciones disponibles"
)
