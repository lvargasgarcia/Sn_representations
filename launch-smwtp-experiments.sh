#!/bin/bash

INSTANCE_DIR="./instances/smtwtp"

filename=$(basename -- "$1")
filename="${filename%.*}"
output="${filename}.esm"

./snob_env/snob_env/bin/python3 fourier.py --mode $2 --problem smwtp --instance "$INSTANCE_DIR/$1" --output $output >&2

cat $output

#Borrar los archivos temporales
rm -rf ./output
