#' @title: Construcción Metadatos
#' @author: Emmanuel Cervantes Monroy
#' @date: Marzo 2026 (Original: 2021)
#' @description: onstruir la mega base de datos uniendo todas las tablas individuales
#   (crudas y derivadas) mediante left_join explícito por px_num.
#   Esta base es el punto de partida para todos los análisis estadísticos
#   posteriores. No se modifica aquí — solo se ensambla y valida.


# Entradas:
#   data/raw/      — bases crudas de 01_importacion.R
#   data/derived/  — variables clínicas de 02_variables_clinicas.R

# Salidas:
#   raw_data.RDS        — Mega base completa
#   colnames.RDA        — Vectores de nombres de columna por tabla

# Anterior: 02_variables_clinicas.R

# Siguiente: cualquier script de análisis (04 en adelante)

# ::::::::::::::::::::::::::::
# Preparación del ambiente
# :::::::::::::::::::::::::::