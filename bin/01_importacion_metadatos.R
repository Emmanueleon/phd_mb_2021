#' @title: Importación Metadatos
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Descargar las bases de datos crudas desde Google Sheets y guardarlas
#  localmente en formato .RDS.


# Entradas:
# 8 hojas de cálculo de Google sheets

# Salidas:
#   mom_demo.RDS     — Demográficos maternos
#   mom_anthro.RDS   — Antropometría materna
#   mom_bx.RDS       — Bioquímicos sanguíneos maternos
#   mom_diet.RDS     — Frecuencia de consumo alimentario materno
#   nb_anthro.RDS    — Antropometría del recién nacido
#   hm_char.RDS      — Características de la leche materna
#   hm_practice.RDS  — Práctica de lactancia
#   feces_char.RDS   — Características de las heces del binomio

# Siguiente Script: 02_variables_clinicas.R


# ::::::::::::::::::::::::::::
# Preparación del ambiente
# :::::::::::::::::::::::::::

## Limpiar el ambiente
rm(list = ls())

## Cargar librerías
suppressPackageStartupMessages({
    library(tidyverse)
    library(gsheet)
    library(here)
})

## Crear carpeta de salida
dir.create(here("data", "raw"), recursive = TRUE, showWarnings = FALSE)

## Definir rutas de salida

out <- list(
    mom_demo    = here("data", "raw", "mom_demo.RDS"),
    mom_anthro  = here("data", "raw", "mom_anthro.RDS"),
    mom_bx      = here("data", "raw", "mom_bx.RDS"),
    mom_diet    = here("data", "raw", "mom_diet.RDS"),
    nb_anthro   = here("data", "raw", "nb_anthro.RDS"),
    hm_char     = here("data", "raw", "hm_char.RDS"),
    hm_practice = here("data", "raw", "hm_practice.RDS"),
    feces_char  = here("data", "raw", "feces_char.RDS")
    )

# :::::::::::::::::::::::::::::
# Cargar datos desde Google Sheets
# ::::::::::::::::::::::::::::::
url <- "https://docs.google.com/spreadsheets/d/1S8f9_z5qfNFXUXCD99AZMrp4Y0spwjQ7PCcG62WFVFQ/edit#gid="

mom_demo <- gsheet2tbl(paste0(url, "0"))
mom_anthro <- gsheet2tbl(paste0(url, "2046898387"))
mom_bx <- gsheet2tbl(paste0(url, "1714043122"))
mom_diet <- gsheet2tbl(paste0(url, "63259731"))
nb_anthro <- gsheet2tbl(paste0(url, "1645407219"))
hm_char <- gsheet2tbl(paste0(url, "1674826407"))
hm_practice <- gsheet2tbl(paste0(url, "1993250597"))
feces_char <- gsheet2tbl(paste0(url, "2018139342"))

# :::::::::::::::::::::::::::::::::::::::::::
# Demográficos de la madre
# ::::::::::::::::::::::::::::::::::::::::::
mom_demo <- mom_demo %>%
    mutate(
        across(c(
            proy, status, clinic, assigned, nacionality, marital.status,
            scolarity, job, pregnancies, vd, csd, abortions, alcohol,
            smoke, drugs, per.history, inhe.history, mom.meds, mom.suppl
        ), as.factor),
        px_num = as.character(px_num),
        no = as.integer(no)
    )

message("mom_demo: ", nrow(mom_demo), " filas, ", ncol(mom_demo), " columnas")
saveRDS(mom_demo, out$mom_demo)

# :::::::::::::::::::::::::::::::::::::::::::
# Antropométricos de la madre
# ::::::::::::::::::::::::::::::::::::::::::
mom_anthro <- mom_anthro %>%
    mutate(across(c(nutc.pre, nutc.imc.2, nutc.fat.2), as.factor))

message("mom_anthro: ", nrow(mom_anthro), " filas, ", ncol(mom_anthro), " columnas")
saveRDS(mom_anthro, out$mom_anthro)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Datos bioquímicos sanguíneos madres
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

message("mom_bx: ", nrow(mom_bx), " filas, ", ncol(mom_bx), " columnas")
saveRDS(mom_bx, out$mom_bx)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Datos dieta materna
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Recodificación de frecuencia de consumo a porciones/día y cálculo de Harris-Benedict
mom_diet <- mom_diet %>%
    mutate(across(
        acelgas:zanahoria,
        ~ case_when(
            .x == 0 ~ 0,
            .x == 1 ~ 0,
            .x == 2 ~ 0.033,
            .x == 3 ~ 0.083,
            .x == 4 ~ 0.143,
            .x == 5 ~ 0.429,
            .x == 6 ~ 0.786,
            .x == 7 ~ 1.0,
            .x == 8 ~ 2.0,
            .x == 9 ~ 2.0,
            TRUE ~ NA_real_ # captura cualquier valor inesperado
        )
    )) %>%
    select(-starts_with("spc"))

message("mom_diet: ", nrow(mom_diet), " filas, ", ncol(mom_diet), " columnas")
saveRDS(mom_diet, out$mom_diet)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 5. Datos antropométricos recién nacido
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Parseo de los datos: clasificación por puntuaciones Z de la OMS
nb_anthro <- nb_anthro %>%
    mutate(across(c(delivery.mode, gender, nb.meds, nb.suppl), as.factor))

message("nb_anthro: ", nrow(nb_anthro), " filas, ", ncol(nb_anthro), " columnas")
saveRDS(nb_anthro, out$nb_anthro)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 6. Datos características leche materna
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
hm_char <- hm_char %>%
    mutate(across(c(hm.1), as.factor))

message("hm_char: ", nrow(hm_char), " filas, ", ncol(hm_char), " columnas")
saveRDS(hm_char, out$hm_char)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 7. Datos práctica lactancia
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Parseo de los datos
hm_practice <- hm_practice %>%
    mutate(across(
        c(
            pre.care, pre.care.freq, pre.care.place, csd.programmed,
            birth.place, birth.contact, birth.meds, medcare.meds,
            medcare, medcare.place, medcare.freq, actual.med,
            colostrum, first.bf, ex.bf, bf.freqday, bf.freqnight,
            bf.extract, bf.type, hm.container, hm.time, hm.temp,
            bf.responsible, bf.others, bf.by.others, co.sleeping,
            activities, activities.place, absence.time, absence.carer,
            bf.carer, dad.job, p1, p1.sex, p1.job, p2, p2.sex, p2.job,
            p3, p3.sex, p3.job, p4, p4.sex, p4.job, p5, p5.sex, p5.job,
            cooking.water, drinking.water, sanitary.service,
            pets, pets.number, pets.type, pets.place, pets.state
        ),
        as.factor
    )) %>%
    mutate(across(
        c(dad.age, p1.age, p2.age, p3.age, p4.age, p5.age),
        as.integer
    ))

message("hm_practice: ", nrow(hm_practice), " filas, ", ncol(hm_practice), " columnas")
saveRDS(hm_practice, out$hm_practice)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 8. Datos características heces
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

feces_char <- feces_char %>%
    mutate(across(
        c(
            mom.st.cons, bristol, mom.st.color,
            nb.st.cons, amsterdam, nb.st.color
        ),
        as.factor
    ))

message("feces_char: ", nrow(feces_char), " filas, ", ncol(feces_char), " columnas")
saveRDS(feces_char, out$feces_char)

# :::::::::::::::::::::::::::::
# Resumen
# :::::::::::::::::::::::::::::
message("\n--- 01_importacion.R completado ---")
message("Archivos guardados en: ", here("data", "raw"))
message("Siguiente paso: correr 02_variables_clinicas.R")