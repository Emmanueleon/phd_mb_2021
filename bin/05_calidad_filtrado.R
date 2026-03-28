#' @title: Calidad filtrado y normalización
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Filtrar y transformar las lecturas del objeto phyloseq crudo
#'   para su posterior analisis. Se realizan los siguientes pasos:
#'   1. Decontam por lotes (primer lote: A,B,C,D,AX,milk / segundo: AZ,BZ,CZ,DZ)
#'   2. Filtrado de muestras con menos de 5000 lecturas
#'   3. Filtrado de ASVs sin lecturas
#'   4. Filtrado de taxas diferentes a Bacteria y Archaea
#'   5. Remocion de controles
#'   6. Transformaciones: composicional, rarefaccion 5k/10k, CLR


# Entradas:
#   out/phyloseq/phy_raw.RDS
#   out/phyloseq/pools_singles.RData

# Salidas:
#   phy_filtered.RDS      - phyloseq filtrado sin transformacion
#   phy_objects.RData     - todos los objetos transformados
#   phy_decontam.RData    - objetos del proceso decontam


# Anterior: 04_importacion_phyloseq.R

# Siguiente: 06_composicion_core.R

# ::::::::::::::::::::::::::::
# Preparación del ambiente
# :::::::::::::::::::::::::::

## Limpiar el ambiente
rm(list = ls())

suppressPackageStartupMessages({
    library(phyloseq)
    library(decontam)
    library(microbiome)
    library(ape)
    library(tidyverse)
    library(here)
})

# ::::::::::::::::::::::::::::::
# ## Definir rutas de salida
# ::::::::::::::::::::::::::::::
out <- list(
    phy_filtered = here("out", "phyloseq", "phy_filtered.RDS"),
    phy_objects  = here("out", "phyloseq", "phy_objects.RData"),
    phy_decontam = here("out", "phyloseq", "phy_decontam.RData")
)

## Cargar datos
message("Cargando datos...")
phy_raw <- readRDS(here("out", "phyloseq", "phy_raw.RDS"))
load(here("out", "phyloseq", "pools_singles.RData"))

message("phy_raw cargado: ", nsamples(phy_raw), " muestras, ", ntaxa(phy_raw), " ASVs")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1. Decontam por lotes
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Se realiza por lotes porque cada corrida de secuenciacion tiene
# sus propios contaminantes ambientales y de reactivos.
# Primer lote:  PoolA, PoolB, PoolC, PoolD, PoolAX, Poolmilk
# Segundo lote: PoolAZ, PoolBZ, PoolCZ, PoolDZ

message("\n--- Decontam ---")

## Extraer OTUs y taxas por lote
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


### Primer lote
message("  Primer lote (A, B, C, D, AX, milk)...")

meta_first <- data.frame(sample_data(phy_raw)) %>%
    filter(rownames(.) %in% c(sample_names(phy_a), sample_names(phy_b),
                               sample_names(phy_c), sample_names(phy_d),
                               sample_names(phy_ax), sample_names(phy_milk))) %>%
    sample_data()

otu_first <- merge_phyloseq(otu_a, otu_b, otu_c, otu_d, otu_milk, otu_ax)
taxa_first <- merge_phyloseq(taxa_a, taxa_b, taxa_c, taxa_d, taxa_milk, taxa_ax)
phy_first <- merge_phyloseq(otu_first, taxa_first, meta_first)

## Identificar contaminantes en primer lote
sample_data(phy_first)$is.neg <- sample_data(phy_first)$sample == "negative"
contamdf_first <- isContaminant(phy_first,
    method = "prevalence",
    neg = "is.neg", threshold = 0.1
)
message("    Contaminantes encontrados: ", sum(contamdf_first$contaminant))

phy_decontam_first <- prune_taxa(!contamdf_first$contaminant, phy_first)
message("    ASVs removidos: ", ntaxa(phy_first) - ntaxa(phy_decontam_first))

### Segundo lote
message("  Segundo lote (AZ, BZ, CZ, DZ)...")

meta_second <- data.frame(sample_data(phy_raw)) %>%
    filter(rownames(.) %in% c(sample_names(phy_az), sample_names(phy_bz),
                               sample_names(phy_cz), sample_names(phy_dz))) %>%
    sample_data()

otu_second <- merge_phyloseq(otu_az, otu_bz, otu_cz, otu_dz)
taxa_second <- merge_phyloseq(taxa_az, taxa_bz, taxa_cz, taxa_dz)
phy_second <- merge_phyloseq(otu_second, taxa_second, meta_second)

## Identificar contaminantes en segundo lote
sample_data(phy_second)$is.neg <- sample_data(phy_second)$sample == "negative"
contamdf_second <- isContaminant(phy_second,
    method = "prevalence",
    neg = "is.neg", threshold = 0.1
)
message("    Contaminantes encontrados: ", sum(contamdf_second$contaminant))

phy_decontam_second <- prune_taxa(!contamdf_second$contaminant, phy_second)
message("    ASVs removidos: ", ntaxa(phy_second) - ntaxa(phy_decontam_second))

### Unir ambos lotes descontaminados
message("  Uniendo lotes descontaminados...")

otu_decontam <- merge_phyloseq(
    otu_table(phy_decontam_first),
    otu_table(phy_decontam_second)
)
taxa_decontam <- merge_phyloseq(
    tax_table(phy_decontam_first),
    tax_table(phy_decontam_second)
)

meta_decontam <- data.frame(sample_data(phy_raw)) %>%
    filter(rownames(.) %in% c(sample_names(phy_decontam_first),
                               sample_names(phy_decontam_second))) %>%
    sample_data()

## Arbol provisional reproducible
set.seed(2021)
tree_decontam <- rtree(nrow(otu_decontam),
    rooted = TRUE,
    tip.label = rownames(otu_decontam)
)

phy_decontam <- merge_phyloseq(
    otu_decontam, taxa_decontam,
    meta_decontam, tree_decontam
)

message(
    "  phy_decontam: ", nsamples(phy_decontam), " muestras, ",
    ntaxa(phy_decontam), " ASVs"
)

phy_object <- phy_decontam


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Filtrado de muestras con menos de 5000 lecturas
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Filtrado de muestras < 5000 lecturas ---")

phy_filsam <- prune_samples(sample_sums(phy_object) > 5000, phy_object)
message("Muestras removidas: ", nsamples(phy_object) - nsamples(phy_filsam))
message("Muestras restantes: ", nsamples(phy_filsam))
phy_object <- phy_filsam


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3. Filtrado de ASVs sin lecturas
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Filtrado de ASVs sin lecturas ---")

phy_filreadtax <- prune_taxa(taxa_sums(phy_object) > 0, phy_object)
message("ASVs removidos: ", ntaxa(phy_object) - ntaxa(phy_filreadtax))
message("ASVs restantes: ", ntaxa(phy_filreadtax))
phy_object <- phy_filreadtax


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4. Filtrado de taxas diferentes a Bacteria y Archaea
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Filtrado de taxas no bacterianas ---")

filter_tax <- c("d__Eukaryota", "Chloroplast", "Mitochondria")

phy_filtax <- subset_taxa(
    phy_object,
    !Kingdom %in% filter_tax &
        !Phylum %in% filter_tax &
        !Class %in% filter_tax &
        !Order %in% filter_tax &
        !Family %in% filter_tax &
        !Genus %in% filter_tax
)
message("ASVs removidos: ", ntaxa(phy_object) - ntaxa(phy_filtax))
message("ASVs restantes: ", ntaxa(phy_filtax))
phy_object <- phy_filtax


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 5. Remocion de controles
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Remocion de controles ---")

proy_samples <- c("mom", "milk", "newborn")
phy_samples <- subset_samples(phy_object, sample %in% proy_samples)

message("Controles removidos: ", nsamples(phy_object) - nsamples(phy_samples))
message("Muestras biologicas: ", nsamples(phy_samples))
message("  ", paste(table(sample_data(phy_samples)$sample),
    names(table(sample_data(phy_samples)$sample)),
    collapse = " | "
))

phy_filtered <- phy_samples


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 6. Transformaciones
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Transformaciones ---")

phy_object <- phy_filtered

## Composicional
message("  Composicional...")
phy_compositional <- microbiome::transform(phy_object, "compositional")

## Rarefaccion
message("  Rarefaccion a 5000 lecturas...")
set.seed(231189)
phy_rarefied5k <- rarefy_even_depth(phy_object,
    sample.size = 5000,
    rngseed = FALSE, verbose = FALSE
)

message("  Rarefaccion a 10000 lecturas...")
set.seed(231189)
phy_rarefied10k <- rarefy_even_depth(phy_object,
    sample.size = 10000,
    rngseed = FALSE, verbose = FALSE
)

## CLR
message("  CLR...")
phy_clr <- microbiome::transform(phy_object, "clr")

message("\n--- Resumen final ---")
message(
    "phy_filtered   : ", nsamples(phy_filtered), " muestras, ",
    ntaxa(phy_filtered), " ASVs"
)
message("phy_rarefied5k : ", nsamples(phy_rarefied5k), " muestras")
message("phy_rarefied10k: ", nsamples(phy_rarefied10k), " muestras")


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 7. Exportacion
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
message("\n--- Exportacion ---")

saveRDS(phy_filtered, out$phy_filtered)
message("phy_filtered guardado")

save(phy_decontam, phy_filsam, phy_filreadtax, phy_filtax, phy_samples,
    phy_compositional, phy_rarefied5k, phy_rarefied10k, phy_clr,
    file = out$phy_objects
)
message("phy_objects guardado")

save(phy_first, phy_decontam_first, contamdf_first,
    phy_second, phy_decontam_second, contamdf_second,
    file = out$phy_decontam
)
message("phy_decontam guardado")

message("\n--- 05_calidad_filtrado.R completado ---")
message("Siguiente paso: 06_composicion_core.R")
