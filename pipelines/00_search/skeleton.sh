#!/bin/bash

# búsqueda de los archivos que acaban con “.fasta” y en lugar de imprimir los errores los mandamos a null.
# en todo el sistema de archivo de la máquina.
printf "Vamos a buscar los archivos tipo fasta."
fastaList=$( find $HOME -name "*.fasta"  2>/dev/null )

printf "$fastaList"

printf "%s\n" "$fastaList" > /tmp/lista-$USER.txt

lineCount=$(printf "%s\n" "$fastaList" | wc -l)

printf "Se encontraron %d archivos fasta.\n" "lineCount"