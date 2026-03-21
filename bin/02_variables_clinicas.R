#' @title: Variable clinicas Metadatos
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: Crear todas las variables clínicas derivadas a partir de los datos crudos
#   generados por 01_importacion.R. Incluye clasificaciones de estado nutricio,
#   puntuaciones Z de la OMS, diagnósticos bioquímicos, ganancia de peso IOM
#   y cálculo del gasto energético basal (Harris-Benedict).

# Entradas:
# Desde data/raw
# mom_anthro.RDS, mom_bx.RDS, mom_diet.RDS, nb_anthro.RDS


# Salidas:
#   mom_anthro_d.RDS   — Antropometría materna + variables derivadas
#   mom_bx_d.RDS       — Bioquímicos + diagnósticos
#   mom_diet_d.RDS     — Dieta + Harris-Benedict
#   nb_anthro_d.RDS    — Antropometría neonatal + puntuaciones Z

# Anterior Script: 01_importacion.R
# Siguiente Script: 03_mega_base.R


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
dir.create(here("data", "derived"), recursive = TRUE, showWarnings = FALSE)

## Cargar datos crudos
mom_anthro <- readRDS(here("data", "raw", "mom_anthro.RDS"))
mom_bx <- readRDS(here("data", "raw", "mom_bx.RDS"))
mom_diet <- readRDS(here("data", "raw", "mom_diet.RDS"))
nb_anthro <- readRDS(here("data", "raw", "nb_anthro.RDS"))


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1. Ingeniería de variables antropometría materna
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
mom_anthro_d <- mom_anthro %>%
    
    ## Ganancia de peso gestacional según criterios IOM
    mutate(iom.gain = case_when(
        nutc.pre == "BP" & between(wt.gain, 12.7, 18.1) ~ "Adecuado",
        nutc.pre == "BP" & wt.gain > 18.1 ~ "Exceso",
        nutc.pre == "BP" & wt.gain < 12.7 ~ "Insuficiente",
        nutc.pre == "NP" & between(wt.gain, 11.3, 15.9) ~ "Adecuado",
        nutc.pre == "NP" & wt.gain > 15.9 ~ "Exceso",
        nutc.pre == "NP" & wt.gain < 11.3 ~ "Insuficiente",
        nutc.pre == "SP" & between(wt.gain, 6.8, 11.3) ~ "Adecuado",
        nutc.pre == "SP" & wt.gain > 11.3 ~ "Exceso",
        nutc.pre == "SP" & wt.gain < 6.8 ~ "Insuficiente",
        nutc.pre == "OB" & between(wt.gain, 5.0, 9.0) ~ "Adecuado",
        nutc.pre == "OB" & wt.gain > 9.0 ~ "Exceso",
        nutc.pre == "OB" & wt.gain < 5.0 ~ "Insuficiente",
        TRUE ~ NA_character_
    )) %>%

    ## Ganancia de peso: 2 categorías
    mutate(iom.gain_2gps = case_when(
        iom.gain == "Adecuado" ~ "Adecuado",
        iom.gain %in% c("Exceso", "Insuficiente") ~ "Inadecuado",
        TRUE ~ NA_character_
    )) %>%

    ## Edad materna punto de corte 30 años
    mutate(mom_age30 = case_when(
        age.calculated <= 30 ~ "young",
        age.calculated > 30 ~ "mature",
        TRUE ~ NA_character_
    )) %>%

    ## Adiposidad punto de corte alternativo 28%
    mutate(nutc.fat28 = case_when(
        bd.fat.2 <= 28 ~ "SN",
        bd.fat.2 > 28 ~ "OB",
        TRUE ~ NA_character_
    )) %>%

    ## Adiposidad punto de corte alternativo 32%
    mutate(nutc.fat32 = case_when(
        bd.fat.2 <= 32 ~ "SN",
        bd.fat.2 > 32 ~ "OB",
        TRUE ~ NA_character_
    )) %>%

    ## Estado nutricio por IMC — 3 grupos
    mutate(nutc_3gps = case_when(
        between(imc.2, 18.5, 24.99) ~ "NP",
        between(imc.2, 25.0, 29.99) ~ "SP",
        imc.2 >= 30 ~ "OB",
        TRUE ~ NA_character_
    )) %>%

    ## Estado nutricio SP+OB combinado (usado en artículo como SYO)
    mutate(nutc.imc.2 = recode(as.character(nutc.imc.2), SP = "SYO", OB = "SYO")) %>%
    
    ## Convertir nuevas variables categóricas a factor
    mutate(across(
        c(
            iom.gain, iom.gain_2gps, mom_age30,
            nutc.fat28, nutc.fat32, nutc_3gps, nutc.imc.2
        ),
        as.factor
    ))

message( "mom_anthro_d: ", ncol(mom_anthro_d) - ncol(mom_anthro), " variables nuevas añadidas")
saveRDS(mom_anthro_d, here("data", "derived", "mom_anthro_d.RDS"))

# -----------------------------------------------------------------------------
# 2. Ingeniería de variables bioquímicos sanguíneos
# -----------------------------------------------------------------------------
mom_bx_d <- mom_bx %>%
    mutate(
        ## Glucosa
        dx.glucose = case_when(
            glucose < 75.0 ~ "Hipoglucemia",
            between(glucose, 75.0, 100.0) ~ "Normoglucemia",
            glucose > 100.0 ~ "Hiperglucemia",
            TRUE ~ NA_character_
        ),

        ## Colesterol total
        dx.cholesterol = case_when(
            cholesterol < 150.0 ~ "Hipocolesterolemia",
            between(cholesterol, 150.0, 200.0) ~ "Normocolesterolemia",
            cholesterol > 200.0 ~ "Hipercolesterolemia",
            TRUE ~ NA_character_
        ),

        ## Triglicéridos
        dx.triglycerides = case_when(
            triglycerides <= 150.0 ~ "Normal",
            triglycerides > 150.0 ~ "Alto",
            TRUE ~ NA_character_
        ),

        ## VLDL
        dx.vldl = case_when(
            vldl < 12.0 ~ "Bajo",
            between(vldl, 12.0, 30.0) ~ "Normal",
            vldl > 30.0 ~ "Alto",
            TRUE ~ NA_character_
        )
    ) %>%
    mutate(across(
        c(dx.glucose, dx.cholesterol, dx.triglycerides, dx.vldl),
        as.factor
    ))

message("mom_bx_d: ", ncol(mom_bx_d) - ncol(mom_bx), " variables nuevas añadidas")
saveRDS(mom_bx_d, here("data", "derived", "mom_bx_d.RDS"))

# -----------------------------------------------------------------------------
# 3. Ingeniería de variables dieta materna Harris-Benedict
# -----------------------------------------------------------------------------
harris_data <- mom_anthro %>%
    select(key, mom.wt.2, mom.ht, mom.age)

mom_diet_d <- mom_diet %>%
    left_join(harris_data, by = "key") %>%
    mutate(Harris = round(
        (10 * mom.wt.2) + (6.25 * mom.ht) - (5 * mom.age) - 161, 2
    )) %>%
    select(-mom.wt.2, -mom.ht, -mom.age) # eliminar columnas auxiliares del join

message(
    "mom_diet_d: Harris-Benedict calculado para ",
    sum(!is.na(mom_diet_d$Harris)), " participantes"
)
saveRDS(mom_diet_d, here("data", "derived", "mom_diet_d.RDS"))

# -----------------------------------------------------------------------------
# 4. Ingeniería de variables antropometría neonatal (puntuaciones Z OMS)
# -----------------------------------------------------------------------------
nb_anthro_d <- nb_anthro %>%
    mutate(
        ## Peso/edad al nacimiento
        z.birth.wa = case_when(
            birth.wa <= -2.0 ~ "Desnutrición",
            birth.wa <= -1.0 & birth.wa > -2.0 ~ "Riesgo desnutrición",
            birth.wa > -1.0 & birth.wa < 1.0 ~ "Eutrófico",
            birth.wa >= 1.0 & birth.wa < 2.0 ~ "Riesgo sobrepeso",
            birth.wa >= 2.0 ~ "Sobrepeso",
            TRUE ~ NA_character_
        ),

        ## Peso/talla al nacimiento
        z.birth.wh = case_when(
            birth.wh < 1.0 ~ "Normal",
            birth.wh >= 1.0 & birth.wh < 2.0 ~ "Sobrepeso",
            birth.wh >= 2.0 ~ "Obesidad",
            TRUE ~ NA_character_
        ),

        ## Talla/edad al nacimiento
        z.birth.ha = case_when(
            birth.ha <= -2.0 ~ "Talla baja",
            birth.ha <= -1.0 & birth.ha > -2.0 ~ "Talla normal baja",
            birth.ha > -1.0 & birth.ha < 1.0 ~ "Normal",
            birth.ha >= 1.0 & birth.ha < 2.0 ~ "Talla normal alta",
            birth.ha >= 2.0 ~ "Talla alta",
            TRUE ~ NA_character_
        ),

        ## Peso/edad al mes
        z.nb.wa = case_when(
            nb.wa <= -2.0 ~ "Desnutrición",
            nb.wa <= -1.0 & nb.wa > -2.0 ~ "Riesgo desnutrición",
            nb.wa > -1.0 & nb.wa < 1.0 ~ "Eutrófico",
            nb.wa >= 1.0 & nb.wa < 2.0 ~ "Riesgo sobrepeso",
            nb.wa >= 2.0 ~ "Sobrepeso",
            TRUE ~ NA_character_
        ),

        ## Peso/talla al mes
        z.nb.wh = case_when(
            nb.wh < 1.0 ~ "Normal",
            nb.wh >= 1.0 & nb.wh < 2.0 ~ "Sobrepeso",
            nb.wh >= 2.0 ~ "Obesidad",
            TRUE ~ NA_character_
        ),

        ## Talla/edad al mes
        z.nb.ha = case_when(
            nb.ha <= -2.0 ~ "Talla baja",
            nb.ha <= -1.0 & nb.ha > -2.0 ~ "Talla normal baja",
            nb.ha > -1.0 & nb.ha < 1.0 ~ "Normal",
            nb.ha >= 1.0 & nb.ha < 2.0 ~ "Talla normal alta",
            nb.ha >= 2.0 ~ "Talla alta",
            TRUE ~ NA_character_
        ),

        ## Perímetro cefálico/edad al mes
        z.nb.cca = case_when(
            nb.cca <= -2.0 ~ "Sospecha microcefalia",
            nb.cca > -2.0 & nb.cca < 2.0 ~ "Normal",
            nb.cca >= 2.0 ~ "Sospecha macrocefalia",
            TRUE ~ NA_character_
        )
    ) %>%
    mutate(across(
        c(
            z.birth.wa, z.birth.wh, z.birth.ha,
            z.nb.wa, z.nb.wh, z.nb.ha, z.nb.cca
        ),
        as.factor
    ))

message(
    "nb_anthro_d: ", ncol(nb_anthro_d) - ncol(nb_anthro),
    " variables nuevas añadidas"
)
saveRDS(nb_anthro_d, here("data", "derived", "nb_anthro_d.RDS"))


# -----------------------------------------------------------------------------
# Resumen final
# -----------------------------------------------------------------------------
message("\n--- 02_variables_clinicas.R completado ---")
message("Archivos guardados en: ", here("data", "derived"))
message("Siguiente paso: correr 03_mega_base.R")



nb_anthro_nb %>% 
select(z.nb.cca) %>% 
print()
