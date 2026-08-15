#!/bin/bash
# ================================================================
# run_weno_tgv.sh — WENO Validation: Taylor-Green Vortex with Shocks
# ================================================================
#
# PURPOSE:
#   Test WENO shock-capturing on the compressible TGV at increasing
#   Mach numbers. At M=0.1 there are no shocks (baseline). At M=0.5
#   weak shocklets form during turbulent transition. At M=1.25
#   strong shocklets appear — C6 should crash, WENO should survive.
#
#   This is the key result: WENO handles shock-turbulence interaction
#   where central schemes fail.
#
# OUTPUTS:
#   weno_tgv_runs/
#   ├── C6_M010/
#   │   ├── input.dat
#   │   ├── log.txt         (full stdout: time, res(1-5), tke, enstrophy)
#   │   └── Monitor.out     (monitor file: time, tke, enstrophy)
#   ├── WENO_JS_M010/
#   ├── ... (one folder per scheme × Mach combo)
#   ├── WENO_CU6M_M125/
#   └── summary.txt
#
# HOW TO RUN:
#   chmod +x run_weno_tgv.sh
#   ./run_weno_tgv.sh
#
# Then open plot_weno_tgv.ipynb in Jupyter.
# ================================================================

set -e

# ================================================================
#  CONFIGURATION
# ================================================================

DSCHEMES=(4      5         6        7          8)
SNAMES=("C6"  "WENO_JS" "WENO_Z" "WENO_CU6" "WENO_CU6M")

# Mach numbers: 0.1 (no shocks), 0.5 (weak shocklets), 1.25 (strong shocklets)
MACHS=(0.1    0.5    1.25)
MLABELS=("M010" "M050" "M125")

# Grid (3D!)
GRID=64
# Uncomment for a quick test (~10× faster, under-resolved but shows trends):
# GRID=32

# Physics
TESTCASE=1        # TGV
VISCOUS=1         # viscous (required for TGV)
RE=1600
GAMMA=1.4
PRANDTL=0.71
T_REF=300

# Time stepping (matching the AS6041 report: dt=0.002, 10000 steps → t=20)
if [ "$GRID" -le 32 ]; then
    DT=0.005;  NSTEPS=4000    # t=20, CFL safe for 32³
else
    DT=0.002;  NSTEPS=10000   # t=20, CFL safe for 64³
fi
RK_STEPS=4

# Filter (F10 with alpha_f=0.495, matching the report)
FSCHEME=10
ALPHA_F=0.495

# Output
RUNDIR="weno_tgv_runs"
ANIMFREQ=99999    # no animation dumps (save disk)

# ================================================================
#  TERMINAL COLOURS
# ================================================================
RED='\033[0;31m'
GRN='\033[0;32m'
YLW='\033[1;33m'
CYN='\033[0;36m'
BLD='\033[1m'
DIM='\033[2m'
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
#  ASCII SPARKLINE for enstrophy after each run
# ================================================================
show_enstrophy_sparkline() {
    local monfile=$1
    if [ ! -f "$monfile" ]; then
        echo -e "  ${DIM}(no monitor file)${RST}"
        return
    fi
    python3 -c "
import sys
chars = '▁▂▃▄▅▆▇█'
try:
    vals = []
    with open('$monfile') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                vals.append(float(parts[2]))  # enstrophy column
    if not vals:
        print('  (empty monitor file)')
        sys.exit(0)
    mn, mx = min(vals), max(vals)
    rng = mx - mn if mx > mn else 1.0
    # Subsample to ~50 chars wide
    step = max(1, len(vals) // 50)
    sampled = vals[::step]
    spark = ''
    for v in sampled:
        idx = int((v - mn) / rng * 7.999)
        idx = max(0, min(7, idx))
        spark += chars[idx]
    peak_val = mx
    peak_idx = vals.index(mx)
    peak_time = peak_idx * 10 * $DT  # every 10 steps
    print(f'  Enstrophy: {spark}  peak={peak_val:.2f} at t≈{peak_time:.1f}')
except Exception as e:
    print(f'  (sparkline error: {e})')
" 2>/dev/null || echo -e "  ${DIM}(sparkline unavailable)${RST}"
}

# ================================================================
#  RESIDUAL TREND after each run
# ================================================================
show_residual_trend() {
    local logfile=$1
    python3 -c "
try:
    vals = []
    with open('$logfile') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 6:
                try:
                    t = float(parts[0])
                    r = max(abs(float(parts[k])) for k in range(1,6))
                    vals.append((t, r))
                except (ValueError, IndexError):
                    continue
    if len(vals) < 2:
        print('  (not enough residual data)')
    else:
        # Check for NaN/Inf
        import math
        nan_found = any(math.isnan(r) or math.isinf(r) for _, r in vals)
        if nan_found:
            print('  Residual: DIVERGED (NaN/Inf detected)')
        else:
            first_r = vals[min(10, len(vals)-1)][1]
            last_r = vals[-1][1]
            if last_r < first_r * 0.1:
                trend = '↓ converging'
            elif last_r > first_r * 10:
                trend = '↑ GROWING'
            else:
                trend = '→ stable'
            print(f'  Residual: {last_r:.2e} ({trend})')
except Exception as e:
    print(f'  (residual parse error: {e})')
" 2>/dev/null || echo -e "  ${DIM}(residual check unavailable)${RST}"
}

# ================================================================
#  COMPILE
# ================================================================
echo ""
echo -e "${BLD}╔══════════════════════════════════════════════════╗${RST}"
echo -e "${BLD}║  WENO TGV VALIDATION — Shocklet Stress Test      ║${RST}"
echo -e "${BLD}╚══════════════════════════════════════════════════╝${RST}"
echo ""
echo -e "${CYN}Compiling solver...${RST}"

gfortran -o acfd -Ofast -march=native -funroll-loops -flto \
  -fdefault-real-8 -fdefault-double-8 \
  Module.f90 compsquare_basic.f90 Postprocessing_routines.f90 \
  Preprocessing_routines.f90 Solver_routines.f90

echo -e "${GRN}✓ Compilation successful.${RST}"
echo ""

# ================================================================
#  SETUP
# ================================================================
TOTAL=$((${#DSCHEMES[@]} * ${#MACHS[@]}))
mkdir -p "$RUNDIR"
COUNT=0
T_START=$(date +%s)

PHYS_TIME=$(echo "$NSTEPS * $DT" | bc)
echo -e "${BLD}Test matrix:${RST} ${#SNAMES[@]} schemes × ${#MACHS[@]} Mach numbers = $TOTAL runs"
echo -e "${BLD}Grid:${RST}        ${GRID}³  ($(echo "$GRID * $GRID * $GRID" | bc) points)"
echo -e "${BLD}Mach:${RST}        ${MACHS[*]}"
echo -e "${BLD}Schemes:${RST}     ${SNAMES[*]}"
echo -e "${BLD}Re:${RST}          $RE    dt=$DT    nsteps=$NSTEPS    t_final=$PHYS_TIME"
echo -e "${BLD}Filter:${RST}      F${FSCHEME}, α_f=$ALPHA_F"
echo ""
echo -e "${YLW}NOTE: 3D TGV at ${GRID}³ is expensive. Estimated ~15-60 min per run.${RST}"
echo -e "${YLW}      Total wall time: ~${TOTAL}×30min ≈ $((TOTAL/2)) hours (rough estimate).${RST}"
echo -e "${YLW}      C6 at M=1.25 is expected to crash — that's the point.${RST}"
echo ""

# ================================================================
#  RUN ALL CASES
# ================================================================

for m in $(seq 0 $((${#MACHS[@]}-1))); do
    mach=${MACHS[$m]}
    mlabel=${MLABELS[$m]}

    for s in $(seq 0 $((${#DSCHEMES[@]}-1))); do
        ds=${DSCHEMES[$s]}
        sname=${SNAMES[$s]}
        COUNT=$((COUNT + 1))

        CASEDIR="${RUNDIR}/${sname}_${mlabel}"
        mkdir -p "$CASEDIR"

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

        echo -e "${BLD}────────────────────────────────────────────────────${RST}"
        printf "  $(progress_bar $((COUNT-1)) $TOTAL)  %s\n" "$eta_str"
        echo -e "  ${YLW}▶ ${sname}${RST} at M=${mach} on ${GRID}³  [$(date '+%H:%M:%S')]"

        # Write input.dat
        cat > input.dat << EOF
0 1
$GRID $GRID $GRID
$TESTCASE $VISCOUS
$RE $mach $GAMMA $PRANDTL $T_REF
6 5
$ds $FSCHEME $ALPHA_F
$RK_STEPS $NSTEPS $DT $ANIMFREQ
EOF
        cp input.dat "$CASEDIR/input.dat"

        # Run solver (don't exit on failure — C6 at high Mach WILL crash)
        set +e
        ./acfd > "$CASEDIR/log.txt" 2>&1
        EXIT_CODE=$?
        set -e

        # Move monitor file
        for f in Monitor_WENO_JS.out Monitor_WENO_Z.out Monitor_WENO_CU6.out \
                 Monitor_WENO_CU6M.out Monitor_C6.out Monitor_E2.out \
                 Monitor_E4.out Monitor_C4.out; do
            [ -f "$f" ] && mv -f "$f" "$CASEDIR/Monitor.out" 2>/dev/null
        done

        # Clean up flow files
        rm -f flow.xyz grid.xyz flow*.xyz

        # Status
        if [ $EXIT_CODE -eq 0 ]; then
            # Check for NaN in the log
            if grep -q "NaN\|Infinity\|nan\|inf" "$CASEDIR/log.txt" 2>/dev/null; then
                echo -e "  ${RED}✗ DIVERGED (NaN detected)${RST}"
            else
                echo -e "  ${GRN}✓ Completed${RST}"
            fi
        else
            echo -e "  ${RED}✗ CRASHED (exit code $EXIT_CODE)${RST}"
        fi

        # Quick diagnostics
        show_residual_trend "$CASEDIR/log.txt"
        show_enstrophy_sparkline "$CASEDIR/Monitor.out"

        echo -e "  ${DIM}Finished: $(date '+%H:%M:%S')  Elapsed: $(elapsed_str $(($(date +%s) - t_now)))${RST}"
    done
done

# ================================================================
#  SUMMARY TABLE
# ================================================================
t_end=$(date +%s)
t_total=$((t_end - T_START))

echo ""
echo -e "${BLD}╔══════════════════════════════════════════════════════════╗${RST}"
echo -e "${BLD}║           WENO TGV — RESULTS SUMMARY                    ║${RST}"
echo -e "${BLD}╚══════════════════════════════════════════════════════════╝${RST}"
echo ""
echo -e "  Total time: ${BLD}$(elapsed_str $t_total)${RST}    $(date)"
echo ""

SUMFILE="${RUNDIR}/summary.txt"
echo "WENO TGV Results — $(date)" > "$SUMFILE"
echo "Grid: ${GRID}³, Re=$RE, Filter: F${FSCHEME} α_f=$ALPHA_F" >> "$SUMFILE"
echo "" >> "$SUMFILE"

# Status matrix
hdr=$(printf "  %-12s" "")
for sname in "${SNAMES[@]}"; do
    hdr="$hdr$(printf '  %-12s' "$sname")"
done
echo "$hdr"
echo "$hdr" >> "$SUMFILE"
div="  $(printf '%-12s' '────────────')"
for sname in "${SNAMES[@]}"; do
    div="$div$(printf '  %-12s' '────────────')"
done
echo "$div"
echo "$div" >> "$SUMFILE"

for m in $(seq 0 $((${#MACHS[@]}-1))); do
    mach=${MACHS[$m]}
    mlabel=${MLABELS[$m]}
    row=$(printf "  M=%-9s" "$mach")
    for sname in "${SNAMES[@]}"; do
        logf="${RUNDIR}/${sname}_${mlabel}/log.txt"
        if [ ! -f "$logf" ]; then
            status="MISSING"
        elif grep -q "NaN\|nan\|Infinity\|inf" "$logf" 2>/dev/null; then
            status="DIVERGED"
        elif grep -q "Writing output done" "$logf" 2>/dev/null; then
            status="✓ OK"
        else
            status="CRASHED"
        fi
        row="$row$(printf '  %-12s' "$status")"
    done
    echo "$row"
    echo "$row" >> "$SUMFILE"
done

echo ""

# Enstrophy peaks
echo "  Enstrophy peaks:"
echo "" >> "$SUMFILE"
echo "  Enstrophy peaks:" >> "$SUMFILE"
hdr2=$(printf "  %-12s" "")
for sname in "${SNAMES[@]}"; do
    hdr2="$hdr2$(printf '  %-12s' "$sname")"
done
echo "$hdr2"
echo "$hdr2" >> "$SUMFILE"
echo "$div"
echo "$div" >> "$SUMFILE"

for m in $(seq 0 $((${#MACHS[@]}-1))); do
    mach=${MACHS[$m]}
    mlabel=${MLABELS[$m]}
    row=$(printf "  M=%-9s" "$mach")
    for sname in "${SNAMES[@]}"; do
        monf="${RUNDIR}/${sname}_${mlabel}/Monitor.out"
        peak=$(python3 -c "
try:
    vals = []
    with open('$monf') as f:
        for line in f:
            p = line.split()
            if len(p) >= 3:
                vals.append(float(p[2]))
    print(f'{max(vals):.2f}' if vals else '—')
except:
    print('—')
" 2>/dev/null)
        row="$row$(printf '  %-12s' "${peak:-—}")"
    done
    echo "$row"
    echo "$row" >> "$SUMFILE"
done

echo ""
echo -e "  ${BLD}What to look for:${RST}"
echo -e "  • M=0.1: all schemes should survive. Compare enstrophy peaks."
echo -e "  • M=0.5: all should survive. WENO peaks may be slightly lower (more dissipation)."
echo -e "  • M=1.25: C6 should ${RED}CRASH${RST}. WENO schemes should ${GRN}SURVIVE${RST}."
echo -e "  •   → This proves WENO handles shock-turbulence interaction."
echo ""
echo -e "  Summary saved to: ${BLD}${RUNDIR}/summary.txt${RST}"
echo -e "  Next: open ${BLD}plot_weno_tgv.ipynb${RST} and run all cells."
echo ""
