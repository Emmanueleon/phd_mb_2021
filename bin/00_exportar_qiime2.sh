#!/usr/bin/env bash
# =============================================================================
# 00_exportar_qiime2.sh
# Proyecto: Microbiota trinomio madre-neonato-leche
# Autor: Emmanuel Cervantes Monroy
# Fecha: Marzo 2026
# =============================================================================
# Objetivo:
#   Debido a que no fue posible usar qiime2R en el ambiente emulado del 2021, se creo este scrip para exportar los artefactos .qza de QIIME2 
#   a formatos legibles por R (.biom, .tsv, .nwk).
#
# Entradas:
#   out/qiime2/Pool*/   — artefactos .qza validados (table, taxonomy, rooted-tree)
#
# Salidas:
#   out/export/Pool*/
#     feature-table.biom  — tabla de ASVs/OTUs
#     taxonomy.tsv        — clasificación taxonómica (SILVA 138)
#     tree.nwk            — árbol filogenético enraizado
#
# Requisitos:
#   conda activate qiime2_2024
#
# Uso:
#   bash bin/00_export_qiime2.sh
#
# NOTAS IMPORTANTES:
#   1. Los .qza fueron generados con QIIME2 2020.6 y validados con QIIME2 2024.2
#      — todos devuelven "valid at level=max"
#   2. PoolD fue secuenciado con un protocolo diferente (Illumina demultiplexado)
#      y no tiene barcodes.fastq.gz — sus .qza son válidos y se usan directamente
#   3. out/export/ está en .gitignore — regenerar corriendo este script
#
#
# Anterior: ninguno (primer script del pipeline de microbiota)
# Siguiente: 04_importacion_phyloseq.R
# =============================================================================

set -e  # detener si cualquier comando falla
set -u  # detener si se usa variable no definida

# -----------------------------------------------------------------------------
# Verificar que QIIME2 está activo
# -----------------------------------------------------------------------------
if ! command -v qiime &> /dev/null; then
    echo "ERROR: QIIME2 no está activo."
    echo "Corre: conda activate qiime2_2024"
    exit 1
fi

echo "QIIME2 version: $(qiime --version)"
echo ""

# -----------------------------------------------------------------------------
# Verificar que estamos en la raíz del proyecto
# -----------------------------------------------------------------------------
if [ ! -d "out/qiime2" ]; then
    echo "ERROR: No se encuentra out/qiime2/"
    echo "Asegúrate de correr este script desde la raíz del proyecto:"
    echo "  cd /Users/manu/Dev/github/Phd_mb_2021"
    exit 1
fi

# -----------------------------------------------------------------------------
# Crear directorios de exportación
# -----------------------------------------------------------------------------
pools=(PoolA PoolB PoolC PoolD PoolAX PoolAZ PoolBZ PoolCZ PoolDZ Poolmilk)

echo "Creando directorios de exportación..."
for pool in "${pools[@]}"; do
    mkdir -p "out/export/$pool"
done
echo "OK"
echo ""

# -----------------------------------------------------------------------------
# Validar todos los .qza antes de exportar
# -----------------------------------------------------------------------------
echo "Validando artefactos .qza..."
validation_errors=0

for qza in out/qiime2/**/*.qza; do
    result=$(qiime tools validate "$qza" 2>&1)
    if echo "$result" | grep -q "valid at level=max"; then
        echo "  OK: $qza"
    else
        echo "  ERROR: $qza"
        echo "  $result"
        validation_errors=$((validation_errors + 1))
    fi
done

if [ $validation_errors -gt 0 ]; then
    echo ""
    echo "ERROR: $validation_errors archivo(s) no válidos. Abortando exportación."
    exit 1
fi

echo ""
echo "Todos los artefactos son válidos. Procediendo con exportación..."
echo ""

# -----------------------------------------------------------------------------
# Exportar cada pool
# -----------------------------------------------------------------------------
for pool in "${pools[@]}"; do

    echo "Exportando $pool..."

    # Obtener sufijo del pool desde el nombre del archivo
    suffix=$(ls out/qiime2/$pool/table*.qza | \
             sed "s|out/qiime2/$pool/table||" | \
             sed 's/\.qza//')

    # Feature table → .biom
    qiime tools export \
        --input-path  "out/qiime2/$pool/table${suffix}.qza" \
        --output-path "out/export/$pool"

    # Taxonomía → taxonomy.tsv
    qiime tools export \
        --input-path  "out/qiime2/$pool/taxonomy${suffix}.qza" \
        --output-path "out/export/$pool"

    # Árbol filogenético → tree.nwk
    qiime tools export \
        --input-path  "out/qiime2/$pool/rooted-tree${suffix}.qza" \
        --output-path "out/export/$pool"

    # Verificar que los 3 archivos se generaron
    file_count=$(ls "out/export/$pool/" | wc -l | tr -d ' ')
    if [ "$file_count" -eq 3 ]; then
        echo "  OK: $file_count archivos generados"
    else
        echo "  ADVERTENCIA: se esperaban 3 archivos, se generaron $file_count"
    fi

done

# -----------------------------------------------------------------------------
# Resumen final
# -----------------------------------------------------------------------------
echo ""
echo "=== Resumen de exportación ==="
total=0
for pool in "${pools[@]}"; do
    count=$(ls "out/export/$pool/" | wc -l | tr -d ' ')
    echo "  $pool: $count archivos"
    total=$((total + count))
done

echo ""
echo "Total: $total archivos exportados"
echo "Directorio: out/export/"
echo ""
echo "--- 00_export_qiime2.sh completado ---"
echo "Siguiente paso: 04_importacion_phyloseq.R"