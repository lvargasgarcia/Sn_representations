#!/bin/bash

INSTANCE_DIR="./instances/smtwtp"
OUTDIR="./output/smwtp/ft/YSR"

mkdir -p "$OUTDIR"  # Crea el directorio de salida si no existe

count=0  # Inicializa el contador

for instance in $INSTANCE_DIR/n5*; do
    filename=$(basename -- "$instance")
    filename="${filename%.*}"

    # Redirige la salida estándar y los errores a /dev/null
    python3 fft.py --problem smwtp --instance "$instance" --output "$OUTDIR/${filename}.ft" > /dev/null 2>&1

    ((count++))  # Incrementa el contador
    if [[ $count -ge 10 ]]; then
        break  # Sale del bucle después de 10 iteraciones
    fi
done

tar -czf - -C "$OUTDIR" .

#Borrar los archivos temporales
rm -rf ./output

