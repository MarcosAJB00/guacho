#!/bin/bash

# entrar a carpeta plots
cd plots || exit

echo "=================================="
echo "Running plotting scripts..."
echo "=================================="

# scripts simples
python inicial.py
python perfil_radial_final.py
python phis_taus.py

# lista de inputs a usar
INPUTS=(1 3 5 8 10)

# loop
for i in "${INPUTS[@]}"
do
    echo "Running evolucion_temporal.py with input $i"
    python evolucion_temporal.py $i
done

echo "=================================="
echo "Plots finished"
echo "=================================="