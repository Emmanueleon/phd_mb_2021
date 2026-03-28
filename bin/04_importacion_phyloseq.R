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
source(here("bin", "funciones_microbiota.R"))


# ::::::::::::::::::::::::::::::::::::::::::::::::
# Construir phyloseq por pool
# ::::::::::::::::::::::::::::::::::::::::::::::::

message("Construyendo objetos phyloseq por pool...")


## Primera corrida::::::::::::::::::::::::::::::::::::::::::::
message("Pool A...")
phy_a <- build_phyloseq("PoolA", here("meta", "metadata_PoolA.tsv"))
phy_a <- prune_samples(c("B116", "M116", "M120"), phy_a)
message("    ", nsamples(phy_a), " muestras conservadas")

message("Pool B...")
phy_b <- build_phyloseq("PoolB", here("meta", "metadata_PoolB.tsv"))
phy_b <- prune_samples(c("B124", "B131", "L134", "M123"), phy_b)
message("    ", nsamples(phy_b), " muestras conservadas")

message("Pool C...")
phy_c <- build_phyloseq("PoolC", here("meta", "metadata_PoolC.tsv"))
phy_c <- prune_samples(c("B68", "L140", "L141", "M66"), phy_c)
message("    ", nsamples(phy_c), " muestras conservadas")

message("Pool D...")
phy_d <- build_phyloseq("PoolD", here("meta", "metadata_PoolD.tsv"))
phy_d <- prune_samples(c(
    "B118", "B130", "B87", "C_Heces", "C_Leche", "H2O",
    "M114", "M122", "M130", "M87", "M94", "MOCK"
), phy_d)
message("    ", nsamples(phy_d), " muestras conservadas (incluye controles)")


## Pool estandarizacion ::::::::::::::::::::::::::::::::::::::::::::
message("Pool Milk...")
phy_milk <- build_phyloseq("Poolmilk", here("meta", "metadata_Poolmilk.tsv"))
phy_milk <- prune_samples(c("L1386B"), phy_milk)
message("    ", nsamples(phy_milk), " muestra conservada")


## Julio 2021 ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("Pool AX...")
phy_ax <- build_phyloseq("PoolAX", here("meta", "metadata_PoolAX.tsv"))
message("    ", nsamples(phy_ax), " muestras conservadas")


## Septiembre 2021 :::::::::::::::::::::::::::::::::::::::::::::::::::::
message("Pool AZ...")
phy_az <- build_phyloseq("PoolAZ", here("meta", "metadata_PoolAZ.tsv"))
phy_az <- prune_samples(setdiff(sample_names(phy_az), c("M66")), phy_az)
if ("M85" %in% sample_names(phy_az)) {
    ## Nota sobre el objeto publicado:
    ## M85 corresponde a muestra materna ("mom").
    sample_data(phy_az)$sample[sample_names(phy_az) == "M85"] <- "mom"
}
message("    ", nsamples(phy_az), " muestras conservadas")
message("    M85 reasignada a 'mom' para reproducibilidad con el objeto publicado")

message("Pool BZ...")
phy_bz <- build_phyloseq("PoolBZ", here("meta", "metadata_PoolBZ.tsv"))
phy_bz <- prune_samples(setdiff(
    sample_names(phy_bz),
    c("B116", "B118", "M116", "M87")
), phy_bz)
## Nota sobre el objeto publicado:
## L131 se conserva en el merge final.
message("    ", nsamples(phy_bz), " muestras conservadas")
message("    L131 conservada para reproducibilidad con el objeto publicado")

message("Pool CZ...")
phy_cz <- build_phyloseq("PoolCZ", here("meta", "metadata_PoolCZ.tsv"))
phy_cz <- phy_cz
## Nota sobre el objeto publicado:
## L131 se conserva en el merge final.
message("    ", nsamples(phy_cz), " muestras conservadas")
message("    L131 conservada para reproducibilidad con el objeto publicado")

message("Pool DZ...")
phy_dz <- build_phyloseq("PoolDZ", here("meta", "metadata_PoolDZ.tsv"))
phy_dz <- prune_samples(setdiff(sample_names(phy_dz), c("M120")), phy_dz)
if ("L117" %in% sample_names(phy_dz)) {
    ## Nota sobre el objeto publicado:
    ## L117 corresponde a muestra de leche ("milk").
    sample_data(phy_dz)$sample[sample_names(phy_dz) == "L117"] <- "milk"
}
message("    ", nsamples(phy_dz), " muestras conservadas")
message("    L117 reasignada a 'milk' para reproducibilidad con el objeto publicado")


# ::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Extraer componentes sin arbol
# ::::::::::::::::::::::::::::::::::::::::::::::::
message("\nExtrayendo componentes por pool...")

otu_a <- otu_table(phy_a)
otu_b <- otu_table(phy_b)
otu_c <- otu_table(phy_c)
otu_d <- otu_table(phy_d)
otu_milk <- otu_table(phy_milk)
otu_ax <- otu_table(phy_ax)
otu_az <- otu_table(phy_az)
otu_bz <- otu_table(phy_bz)
otu_cz <- otu_table(phy_cz)
otu_dz <- otu_table(phy_dz)

taxa_a <- tax_table(phy_a)
taxa_b <- tax_table(phy_b)
taxa_c <- tax_table(phy_c)
taxa_d <- tax_table(phy_d)
taxa_milk <- tax_table(phy_milk)
taxa_ax <- tax_table(phy_ax)
taxa_az <- tax_table(phy_az)
taxa_bz <- tax_table(phy_bz)
taxa_cz <- tax_table(phy_cz)
taxa_dz <- tax_table(phy_dz)


# ::::::::::::::::::::::::::::::::::::::::::::::::
# 3. Merge de OTUs y taxonomias
# ::::::::::::::::::::::::::::::::::::::::::::::::
message("Haciendo merge de OTUs y taxonomias...")

otu_raw <- merge_phyloseq(
    otu_a, otu_b, otu_c, otu_d,
    otu_ax, otu_milk,
    otu_az, otu_bz, otu_cz, otu_dz
)

taxa_raw <- merge_phyloseq(
    taxa_a, taxa_b, taxa_c, taxa_d,
    taxa_ax, taxa_milk,
    taxa_az, taxa_bz, taxa_cz, taxa_dz
)

message("ASVs totales: ", nrow(otu_raw))


# ::::::::::::::::::::::::::::::::::::::::::::::::
# 4. Metadatos conjuntos
# ::::::::::::::::::::::::::::::::::::::::::::::::
message("Construyendo metadatos conjuntos...")

# Definir columnas estándar
cols_standard <- c(
    "barcode.sequence", "match", "place", "sample",
    "imc", "fat", "delivery", "sex"
)

# Wrangling: quedarse solo con columnas estándar en cada pool
meta_list <- lapply(list(
    as.data.frame(sample_data(phy_a)),
    as.data.frame(sample_data(phy_b)),
    as.data.frame(sample_data(phy_c)),
    as.data.frame(sample_data(phy_d)),
    as.data.frame(sample_data(phy_ax)),
    as.data.frame(sample_data(phy_milk)),
    as.data.frame(sample_data(phy_az)),
    as.data.frame(sample_data(phy_bz)),
    as.data.frame(sample_data(phy_cz)),
    as.data.frame(sample_data(phy_dz))
), function(x) x[, cols_standard])

# Merge limpio
meta_raw <- do.call(rbind, meta_list)
meta <- sample_data(meta_raw)
message("  Muestras totales: ", nrow(meta_raw))

# ::::::::::::::::::::::::::::::::::::::::::::::::
# 5. Arbol provisional reproducible
# ::::::::::::::::::::::::::::::::::::::::::::::::
message("Generando arbol provisional (set.seed = 2021)...")
message("PENDIENTE: sustituir por arbol real al rehacer QIIME2 desde fastq.gz")

set.seed(2021)
tree_raw <- rtree(
    ntaxa(otu_raw),
    rooted    = TRUE,
    tip.label = taxa_names(otu_raw)
)


# ::::::::::::::::::::::::::::::::::::::::::::::::
# 6. Objeto phyloseq final
# ::::::::::::::::::::::::::::::::::::::::::::::::
message("Ensamblando objeto phyloseq final...")

phy_raw <- merge_phyloseq(otu_raw, taxa_raw, meta, tree_raw)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 7. Enriquecer sample_data con variables derivadas
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("Añadiendo variables derivadas al sample_data...")

## Cargar master_data para join
master_data <- readRDS(here("data", "processed", "master_data.RDS"))

## Seleccionar variables interes
vars_extra <- master_data %>%
    select(px_num, mom_age30, nutc_3gps, iom.gain) %>%
    rename(match = px_num) 


## Crear  data frame del phy_raw
meta_enriquecida <- data.frame(sample_data(phy_raw)) %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    left_join(vars_extra, by = "match") %>%
    tibble::column_to_rownames("sample_id")

## Reasignar el phyloseq al nuevo sample_data
sample_data(phy_raw) <- sample_data(meta_enriquecida)

message("\n--- Validacion phy_raw ---")
message("Muestras : ", nsamples(phy_raw))
message("ASVs     : ", ntaxa(phy_raw))
message("Variables: ", length(sample_variables(phy_raw)))
message("Variables en sample_data: ",
        paste(colnames(sample_data(phy_raw)), collapse = ", "))


# ::::::::::::::::::::::::::::::::::::::::::::::::
# 7. Exportacion
# ::::::::::::::::::::::::::::::::::::::::::::::::
saveRDS(phy_raw, out$phy_raw)
message("\nphy_raw guardado: ", out$phy_raw)

save(phy_a, phy_b, phy_c, phy_d,
    phy_ax, phy_milk,
    phy_az, phy_bz, phy_cz, phy_dz,
    file = out$pools_singles
)
message("pools_singles guardado")

message("\n--- 04_importacion_phyloseq.R completado ---")
message("Siguiente paso: 05_calidad_filtrado.R")
message("\nPENDIENTES:")
message("  [ ] Documentar que L131 se conserva para reproducir el objeto publicado")
message("  [ ] Sustituir arbol provisional al rehacer QIIME2 desde fastq.gz")
