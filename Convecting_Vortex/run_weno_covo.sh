#!/bin/bash
# ================================================================
# run_weno_covo.sh — WENO Validation via COVO Convergence Study
# ================================================================
#
# PURPOSE:
#   Test all 4 WENO variants (JS, Z, CU6, CU6-M) on the convecting
#   vortex (COVO) and compare against the C6 compact baseline.
#   This tells you:
#     - Does each WENO scheme converge at the right order?
#     - Does it match/beat C6 on smooth problems?
#     - Is the vortex shape preserved?
#
# WHAT THIS DOES:
#   1. Compiles the solver
#   2. Runs 15 simulations: 5 schemes × 3 grids (uniform COVO)
#   3. Optionally runs on the sinusoidal mesh for curvilinear check
#   4. Organises output into weno_covo_runs/ with subfolders
#   5. Prints a summary table with errors and convergence rates
#
# FOLDER STRUCTURE:
#   weno_covo_runs/
#   ├── C6_N041/
#   │   ├── input.dat
#   │   ├── log.txt
#   │   ├── covo_centreline.dat
#   │   └── covo_vorticity.dat
#   ├── C6_N081/
#   ├── ...
#   ├── WENO_JS_N041/
#   ├── WENO_JS_N081/
#   ├── ...
#   ├── WENO_CU6M_N161/
#   └── summary.txt
#
# HOW TO RUN:
#   chmod +x run_weno_covo.sh
#   ./run_weno_covo.sh
#
# Then open plot_weno_covo.ipynb in Jupyter and run all cells.
# ================================================================

set -e

# ================================================================
#  CONFIGURATION — edit these if needed
# ================================================================

# Schemes: C6 baseline + 4 WENO variants
DSCHEMES=(4      5         6        7          8)
SNAMES=("C6"  "WENO_JS" "WENO_Z" "WENO_CU6" "WENO_CU6M")

# Grid sizes for convergence study
GRIDS=(41 81 161)
# Uncomment the line below to add 321 (takes much longer):
# GRIDS=(41 81 161 321)

# Physics (inviscid COVO at M=0.1)
TESTCASE=2        # 2 = uniform COVO
VISCOUS=0
RE=1600
MACH=0.1
GAMMA=1.4
PRANDTL=0.71
T_REF=300

# Numerics
DT=0.002          # time step (CFL safe for all grids up to 321)
NSTEPS=8000       # total steps → physical time = 16.0 (one full lap)
RK_STEPS=4

# Filter (Visbal–Gaitonde SF6, mild)
FSCHEME=6
ALPHA_F=0.3

# Output
RUNDIR="weno_covo_runs"
ANIMFREQ=99999    # no animation dumps (saves disk)

# ================================================================
#  CURVILINEAR TEST (optional)
#  Set to 1 to also run on the sinusoidal mesh (testcase 4).
#  This tests that the curvilinear metrics are handled correctly.
#  Only runs on a single grid (81) with all schemes.
# ================================================================
RUN_CURVILINEAR=0
CURV_GRID=81

# ================================================================
#  HELPER FUNCTIONS
# ================================================================

# Colours
RED='\033[0;31m'
GRN='\033[0;32m'
YLW='\033[1;33m'
CYN='\033[0;36m'
BLD='\033[1m'
RST='\033[0m'

elapsed_str() {
    local s=$1
    printf '%02d:%02d:%02d' $((s/3600)) $(( (s%3600)/60 )) $((s%60))
}

progress_bar() {
    local done=$1 total=$2 width=30
    local filled=$((done * width / total))
    local empty=$((width - filled))
    printf '['
    printf '%0.s█' $(seq 1 $filled 2>/dev/null) || true
    printf '%0.s░' $(seq 1 $empty 2>/dev/null)  || true
    printf '] %d/%d' "$done" "$total"
}

# ================================================================
#  STEP 1: COMPILE
# ================================================================
echo ""
echo -e "${BLD}╔══════════════════════════════════════════════╗${RST}"
echo -e "${BLD}║   WENO COVO VALIDATION — CompSquare Solver   ║${RST}"
echo -e "${BLD}╚══════════════════════════════════════════════╝${RST}"
echo ""
echo -e "${CYN}Compiling solver...${RST}"

gfortran -o acfd -Ofast -march=native -funroll-loops -flto \
  -fdefault-real-8 -fdefault-double-8 \
  Module.f90 compsquare_basic.f90 Postprocessing_routines.f90 \
  Preprocessing_routines.f90 Solver_routines.f90

echo -e "${GRN}✓ Compilation successful.${RST}"
echo ""

# ================================================================
#  STEP 2: COUNT TOTAL RUNS AND SET UP
# ================================================================
TOTAL=$((${#DSCHEMES[@]} * ${#GRIDS[@]}))
if [ "$RUN_CURVILINEAR" -eq 1 ]; then
    TOTAL=$((TOTAL + ${#DSCHEMES[@]}))
fi

mkdir -p "$RUNDIR"
COUNT=0
T_START=$(date +%s)

echo -e "${BLD}Test matrix:${RST} ${#DSCHEMES[@]} schemes × ${#GRIDS[@]} grids = $TOTAL runs"
echo -e "${BLD}Schemes:${RST}     ${SNAMES[*]}"
echo -e "${BLD}Grids:${RST}       ${GRIDS[*]}"
echo -e "${BLD}Phys. time:${RST}  $(echo "$NSTEPS * $DT" | bc) (one full vortex lap)"
if [ "$RUN_CURVILINEAR" -eq 1 ]; then
    echo -e "${BLD}Curvilinear:${RST} ON (sinusoidal mesh, N=${CURV_GRID})"
fi
echo ""

# ================================================================
#  STEP 3: RUN ALL UNIFORM COVO CASES
# ================================================================

run_case() {
    local ds=$1 sname=$2 ni=$3 tc=$4 label=$5
    local npad=$(printf '%03d' $ni)
    local casedir="${RUNDIR}/${label}_N${npad}"
    mkdir -p "$casedir"

    # Write input.dat
    cat > input.dat << EOF
0 1
$ni $ni 1
$tc $VISCOUS
$RE $MACH $GAMMA $PRANDTL $T_REF
6 5
$ds $FSCHEME $ALPHA_F
$RK_STEPS $NSTEPS $DT $ANIMFREQ
EOF
    cp input.dat "$casedir/input.dat"

    # Run solver
    ./acfd > "$casedir/log.txt" 2>&1
    local exit_code=$?

    # Move output files
    for f in Monitor_COVO_N${npad}.out Monitor_COVO_RAND_N${npad}.out Monitor_COVO_SINE_N${npad}.out; do
        [ -f "$f" ] && mv -f "$f" "$casedir/Monitor.out" 2>/dev/null
    done
    [ -f "covo_centreline_N${npad}.dat" ] && mv -f "covo_centreline_N${npad}.dat" "$casedir/covo_centreline.dat"
    [ -f "covo_vorticity_N${npad}.dat" ]  && mv -f "covo_vorticity_N${npad}.dat"  "$casedir/covo_vorticity.dat"
    rm -f flow.xyz grid.xyz flow*.xyz

    return $exit_code
}

for s in $(seq 0 $((${#DSCHEMES[@]}-1))); do
    ds=${DSCHEMES[$s]}
    sname=${SNAMES[$s]}

    for ni in "${GRIDS[@]}"; do
        COUNT=$((COUNT + 1))
        npad=$(printf '%03d' $ni)

        # Progress
        t_now=$(date +%s)
        t_elapsed=$((t_now - T_START))
        if [ $COUNT -gt 1 ]; then
            t_per_run=$((t_elapsed / (COUNT - 1)))
            t_remaining=$(( t_per_run * (TOTAL - COUNT + 1) ))
            eta_str="ETA $(elapsed_str $t_remaining)"
        else
            eta_str=""
        fi

        echo -e "${BLD}────────────────────────────────────────────${RST}"
        printf "  $(progress_bar $((COUNT-1)) $TOTAL)  %s\n" "$eta_str"
        echo -e "  ${YLW}▶ ${sname}${RST} on ${ni}×${ni}  [$(date '+%H:%M:%S')]"

        if run_case "$ds" "$sname" "$ni" "$TESTCASE" "$sname"; then
            # Extract error from log
            err=$(grep "COVO L-inf" "${RUNDIR}/${sname}_N${npad}/log.txt" 2>/dev/null | awk '{print $NF}')
            if [ -n "$err" ]; then
                echo -e "  ${GRN}✓${RST} L∞ = ${BLD}${err}${RST}"
            else
                echo -e "  ${GRN}✓${RST} Completed (no error line found)"
            fi
        else
            echo -e "  ${RED}✗ CRASHED — check ${RUNDIR}/${sname}_N${npad}/log.txt${RST}"
        fi
    done
done

# ================================================================
#  STEP 4: OPTIONAL CURVILINEAR TEST
# ================================================================
if [ "$RUN_CURVILINEAR" -eq 1 ]; then
    echo ""
    echo -e "${BLD}──── Curvilinear mesh test (sinusoidal, testcase 4) ────${RST}"
    for s in $(seq 0 $((${#DSCHEMES[@]}-1))); do
        ds=${DSCHEMES[$s]}
        sname=${SNAMES[$s]}
        COUNT=$((COUNT + 1))
        label="${sname}_SINE"
        npad=$(printf '%03d' $CURV_GRID)

        printf "  $(progress_bar $((COUNT-1)) $TOTAL)  "
        echo -e "${YLW}▶ ${sname}${RST} on sinusoidal ${CURV_GRID}×${CURV_GRID}"

        if run_case "$ds" "$sname" "$CURV_GRID" 4 "$label"; then
            err=$(grep "COVO L-inf" "${RUNDIR}/${label}_N${npad}/log.txt" 2>/dev/null | awk '{print $NF}')
            [ -n "$err" ] && echo -e "  ${GRN}✓${RST} L∞ = ${BLD}${err}${RST}" || echo -e "  ${GRN}✓${RST} Done"
        else
            echo -e "  ${RED}✗ CRASHED${RST}"
        fi
    done
fi

# ================================================================
#  STEP 5: SUMMARY TABLE
# ================================================================
t_end=$(date +%s)
t_total=$((t_end - T_START))

echo ""
echo -e "${BLD}╔══════════════════════════════════════════════════════════╗${RST}"
echo -e "${BLD}║              WENO COVO — RESULTS SUMMARY                ║${RST}"
echo -e "${BLD}╚══════════════════════════════════════════════════════════╝${RST}"
echo ""
echo -e "  Total time: ${BLD}$(elapsed_str $t_total)${RST}    $(date)"
echo ""

SUMFILE="${RUNDIR}/summary.txt"
echo "WENO COVO Convergence Results — $(date)" > "$SUMFILE"
echo "" >> "$SUMFILE"

# Header
hdr=$(printf "  %-12s" "Grid")
for sname in "${SNAMES[@]}"; do
    hdr="$hdr$(printf '  %-14s' "$sname")"
done
div="  $(printf '%-12s' '────────────')"
for sname in "${SNAMES[@]}"; do
    div="$div$(printf '  %-14s' '──────────────')"
done

echo "$hdr"
echo "$div"
echo "$hdr" >> "$SUMFILE"
echo "$div" >> "$SUMFILE"

for ni in "${GRIDS[@]}"; do
    npad=$(printf '%03d' $ni)
    row=$(printf "  %-12s" "${ni}×${ni}")
    for sname in "${SNAMES[@]}"; do
        logf="${RUNDIR}/${sname}_N${npad}/log.txt"
        err=$(grep "COVO L-inf" "$logf" 2>/dev/null | awk '{print $NF}')
        if [ -z "$err" ]; then
            row="$row$(printf '  %-14s' '—')"
        else
            row="$row$(printf '  %-14s' "$err")"
        fi
    done
    echo "$row"
    echo "$row" >> "$SUMFILE"
done

# Convergence rates (computed between consecutive grids)
echo ""
echo "  Convergence rates (between consecutive grids):"
echo "" >> "$SUMFILE"
echo "  Convergence rates:" >> "$SUMFILE"

prev_grids=()
for idx in $(seq 1 $((${#GRIDS[@]}-1))); do
    ni_coarse=${GRIDS[$((idx-1))]}
    ni_fine=${GRIDS[$idx]}
    dx_coarse=$(echo "16.0 / ($ni_coarse - 1)" | bc -l)
    dx_fine=$(echo "16.0 / ($ni_fine - 1)" | bc -l)
    ratio_label="${ni_coarse}→${ni_fine}"

    row=$(printf "  %-12s" "$ratio_label")
    for sname in "${SNAMES[@]}"; do
        npad_c=$(printf '%03d' $ni_coarse)
        npad_f=$(printf '%03d' $ni_fine)
        err_c=$(grep "COVO L-inf" "${RUNDIR}/${sname}_N${npad_c}/log.txt" 2>/dev/null | awk '{print $NF}')
        err_f=$(grep "COVO L-inf" "${RUNDIR}/${sname}_N${npad_f}/log.txt" 2>/dev/null | awk '{print $NF}')
        if [ -n "$err_c" ] && [ -n "$err_f" ]; then
            rate=$(python3 -c "
import math
e1, e2 = $err_c, $err_f
d1, d2 = $dx_coarse, $dx_fine
if e1 > 0 and e2 > 0 and d1 > 0 and d2 > 0:
    print(f'{math.log(e1/e2)/math.log(d1/d2):.2f}')
else:
    print('—')
" 2>/dev/null)
            row="$row$(printf '  %-14s' "$rate")"
        else
            row="$row$(printf '  %-14s' '—')"
        fi
    done
    echo "$row"
    echo "$row" >> "$SUMFILE"
done

echo ""
echo -e "  Expected rates:  C6 ≈ 4–6    JS ≈ 5    Z ≈ 5    CU6 ≈ 6    CU6-M ≈ 6"
echo ""

if [ "$RUN_CURVILINEAR" -eq 1 ]; then
    echo "  Curvilinear (sinusoidal mesh, ${CURV_GRID}×${CURV_GRID}):"
    for sname in "${SNAMES[@]}"; do
        npad=$(printf '%03d' $CURV_GRID)
        err=$(grep "COVO L-inf" "${RUNDIR}/${sname}_SINE_N${npad}/log.txt" 2>/dev/null | awk '{print $NF}')
        printf "    %-12s  %s\n" "$sname" "${err:-(not found)}"
    done
    echo ""
fi

echo -e "  Summary saved to: ${BLD}${RUNDIR}/summary.txt${RST}"
echo -e "  Next: open ${BLD}plot_weno_covo.ipynb${RST} and run all cells."
echo ""
