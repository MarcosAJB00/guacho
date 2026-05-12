#!/bin/bash

# entrar a carpeta plots
cd plots || exit

echo "Running plotting scripts..."

# scripts simples
python inicial.py
python perfil_radial_final.py
python phis_taus.py

# lista de inputs a usar
INPUTS=(10 200 500 1000 1250 1500 1750 2000)

# loop
for i in "${INPUTS[@]}"
do
    echo "Running evolucion_temporal.py with input $i"
    python evolucion_temporal.py $i
done

echo "Plots finished"