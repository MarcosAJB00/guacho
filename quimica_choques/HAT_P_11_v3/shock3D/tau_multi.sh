#!/bin/bash

for d in phi*; do
     echo "Procesando $d"
    (
      cd "$d" || exit
      python prueba_transito_2.py
    )
done
