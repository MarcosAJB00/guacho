#!/bin/bash
#
# run_shock3D.sh
#
# Automatiza, para cada carpeta phi<N>_alpha<M>/:
#   1) parchea phi_c y alpha en plot.py con los valores del nombre de carpeta
#   2) corre `python plot.py`
#   3) encola submit_tau0.sh, submit_tau1.sh y submit_tau2.sh con sbatch
#
# Uso:
#   ./run_shock3D.sh            # corre todo en serio
#   ./run_shock3D.sh --dry-run  # solo muestra qué haría, no modifica ni corre nada
#
set -eu
set -o pipefail 2>/dev/null || true   # algunos bash viejos no soportan pipefail, se ignora si falla

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "*** MODO DRY-RUN: no se modifica nada ni se encola nada ***"
fi

BASE_DIR="$(pwd)"

for dir in phi*_alpha*; do
    [[ -d "$dir" ]] || continue

    # Extrae los números de phi y alpha del nombre de la carpeta
    if [[ "$dir" =~ ^phi([0-9]+)_alpha([0-9]+)$ ]]; then
        phi="${BASH_REMATCH[1]}"
        alpha="${BASH_REMATCH[2]}"
    else
        echo ">> Carpeta '$dir' no matchea el patrón phi<N>_alpha<M>, se salta."
        continue
    fi

    echo "=============================================="
    echo ">> Procesando $dir  (phi_c=${phi}.0, alpha=${alpha}.0)"
    echo "=============================================="

    cd "$BASE_DIR/$dir"

    if [[ ! -f plot.py ]]; then
        echo "   !! No se encontró plot.py en $dir, se salta."
        cd "$BASE_DIR"
        continue
    fi

    # Backup del plot.py original (solo la primera vez)
    if [[ ! -f plot.py.orig ]]; then
        cp plot.py plot.py.orig
    fi

    if $DRY_RUN; then
        echo "   [dry-run] sed reemplazaría phi_c -> ${phi}.0"
        echo "   [dry-run] sed reemplazaría alpha -> ${alpha}.0"
        echo "   [dry-run] python plot.py"
        echo "   [dry-run] sbatch submit_tau0.sh submit_tau1.sh submit_tau2.sh"
        cd "$BASE_DIR"
        continue
    fi

    # --- Parchea phi_c ---
    n_phi=$(grep -cE '^phi_c[[:space:]]*=[[:space:]]*np\.deg2rad\(' plot.py || true)
    if [[ "$n_phi" -eq 0 ]]; then
        echo "   !! No se encontró línea de phi_c en $dir/plot.py, reviso manualmente."
        cd "$BASE_DIR"
        continue
    fi
    sed -i -E "s/(^phi_c[[:space:]]*=[[:space:]]*np\.deg2rad\()[0-9.]+(\))/\1${phi}.0\2/" plot.py

    # --- Parchea alpha ---
    n_alpha=$(grep -cE '^alpha[[:space:]]*=[[:space:]]*np\.deg2rad\(' plot.py || true)
    if [[ "$n_alpha" -eq 0 ]]; then
        echo "   !! No se encontró línea de alpha en $dir/plot.py, reviso manualmente."
        cd "$BASE_DIR"
        continue
    fi
    sed -i -E "s/(^alpha[[:space:]]*=[[:space:]]*np\.deg2rad\()[0-9.]+(\))/\1${alpha}.0\2/" plot.py

    echo "   Valores parcheados. Verificación:"
    grep -E '^(phi_c|alpha)[[:space:]]*=' plot.py

    # --- Corre plot.py ---
    echo "   Corriendo python plot.py ..."
    python plot.py

    # --- Encola los jobs ---
    for job in submit_tau0.sh submit_tau1.sh submit_tau2.sh; do
        if [[ -f "$job" ]]; then
            echo "   sbatch $job"
            sbatch "$job"
        else
            echo "   !! No se encontró $job en $dir, se salta ese sbatch."
        fi
    done

    cd "$BASE_DIR"
done

echo "=============================================="
echo "Listo."