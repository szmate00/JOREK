#!/bin/bash
#
# transport_fit.sh -- unattended fitting of the tabulated D_perp / ZK_perp profiles of a JOREK
# model600 run (kinetic_main) so that the upstream profiles hold the step-0 profiles and the
# target profiles follow Langmuir-probe data.  See doc/transport_fit/README.md.
#
# No command-line arguments: every setting is a variable below.  Everything lives flat in $BASE:
#   $BASE/$NAMELIST            the run's namelist (restart, nstep_n, restart_particles are set here)
#   $BASE/$EXE                 the executable (kinetic_main)
#   $BASE/jorek2_postproc      the postprocessor
#   $BASE/transport_fit.py     the python core (copy from util/)
#   $BASE/transport_fit_toy.py only for RUN_MODE=toy
#   $BASE/probes.dat           R Z ne[1e19 m-3] Te[eV] [ne_err Te_err], one probe per line
#   $BASE/jorek??????.h5       step files; the loop restarts from the newest one, REF_STEP is the reference
# The loop writes tfit.cfg, tfit_state.json, tfit_log.csv, tfit_iterNN.png, tfit_driver.log, the
# postproc command files tfit_pp_*.pp and backups <profile file>.iterNN / .user.
# Re-running the script resumes where it stopped (tfit_state.json).
#
# ============================================ settings ============================================
BASE=/tokp/work/szmate/tfit_case      # run directory
NAMELIST=input
EXE=kinetic_main                      # executable in $BASE
RUN_MODE=slurm                        # slurm | local | toy
PROBE_FILE=probes.dat                 # relative to $BASE; empty = upstream-only fit (targets reported)
TEMPERATURES=2                        # 1: T with ZK_perp_file; 2: T_e, T_i with ZK_e_perp_file, ZK_i_perp_file
REF_STEP=0                            # upstream reference = jorek<REF_STEP>.h5 (the initial profiles)
NSTEP_SEG=200                         # steps per segment (doubled when the state is not steady)
NSTEP_MAX=800
MAX_ITER=12
TOL=0.05                              # converged: all channels within TOL (rel. rms) for two iterations
RELAX=0.7                             # 0.5 with drifts
CLAMP=2.0                             # max factor per knot and iteration; 1.5 with drifts
PSI_MIN=0.7                           # knots below are held (the core equilibrates too slowly)
W_TARGET=0.5                          # weight of the target width in the SOL knots (0 = upstream only)
TE_TARGET_WEIGHT=0.5                  # extra factor on W_TARGET for the probe T_e
CHI_EDGE_MAX=10                       # m2/s cap on chi near the grid edge for the first three iterations
LINEPOINTS=2000                       # points of the outer-midplane cut
EXTRA_CFG=""                          # extra KEY=VALUE pairs for tfit.cfg, e.g. "PLATE_ZMAX=-0.5 KNOTS=0,0.5,0.8,..."
# --- slurm (RUN_MODE=slurm)
JOB_NAME=tfit
PARTITION=p.tok
QOS=p.tok.2h
NODES=2
NTASKS_PER_NODE=8
CPUS_PER_TASK=32
WALLTIME=00:20:00
MAIL=mate.szuecs.extern@ipp.mpg.de
# --- local (RUN_MODE=local)
NP=8
# ======================================== end of settings ========================================

cd "$BASE" || { echo "BASE=$BASE does not exist"; exit 1; }
PY=$BASE/transport_fit.py
TOY=$BASE/transport_fit_toy.py
POSTPROC=$BASE/jorek2_postproc
CFG=$BASE/tfit.cfg
LOG=$BASE/tfit_driver.log

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

newest_step() {
  ls jorek[0-9]*.h5 2>/dev/null | sed -E 's/jorek0*([0-9]+)\.h5/\1/' | sort -n | tail -1
}

write_cfg() {
  {
    echo "RUN=$BASE"
    echo "NAMELIST=$NAMELIST"
    echo "TEMPERATURES=$TEMPERATURES"
    echo "PROBE_FILE=$PROBE_FILE"
    echo "REF_STEP=$REF_STEP"
    echo "NSTEP_SEG=$NSTEP_SEG"
    echo "NSTEP_MAX=$NSTEP_MAX"
    echo "MAX_ITER=$MAX_ITER"
    echo "TOL=$TOL"
    echo "RELAX=$RELAX"
    echo "CLAMP=$CLAMP"
    echo "PSI_MIN=$PSI_MIN"
    echo "W_TARGET=$W_TARGET"
    echo "TE_TARGET_WEIGHT=$TE_TARGET_WEIGHT"
    echo "CHI_EDGE_MAX=$CHI_EDGE_MAX"
    echo "LINEPOINTS=$LINEPOINTS"
    for kv in $EXTRA_CFG; do echo "$kv"; done
  } > "$CFG"
}

write_job() {
  cat > tfit_job.sh <<EOF
#!/bin/bash -l
#SBATCH -J $JOB_NAME
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#SBATCH -D ./
#SBATCH --partition=$PARTITION
#SBATCH --qos=$QOS
#SBATCH --nodes=$NODES
#SBATCH --ntasks-per-node=$NTASKS_PER_NODE
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#SBATCH --mem 0
#SBATCH --time=$WALLTIME
#SBATCH --mail-type=fail
#SBATCH --mail-user=$MAIL

ulimit -s unlimited
export OMP_STACKSIZE="500M"
export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK:-1}
export OMP_PLACES=cores
set -o pipefail
srun ./$EXE < $NAMELIST | tee logfile
EOF
}

run_segment() {                      # one restart segment of $NSTEP steps; non-zero = failed
  case "$RUN_MODE" in
    slurm) write_job; sbatch --wait tfit_job.sh ;;
    local) mpirun -np "$NP" "./$EXE" < "$NAMELIST" > logfile 2>&1 ;;
    toy)   python3 "$TOY" run --cfg "$CFG" --nstep "$NSTEP" ;;
    *)     log "unknown RUN_MODE=$RUN_MODE"; exit 1 ;;
  esac
}

postproc_step() {                    # exports for one step: outer-midplane cut and boundary quantities
  local step=$1
  if [ "$RUN_MODE" = toy ]; then
    python3 "$TOY" export --cfg "$CFG" --step "$step"
    return $?
  fi
  cat > "tfit_pp_$step.pp" <<EOF
namelist $NAMELIST
si-units
set linepoints $LINEPOINTS
expressions $MIDPLANE
mark_coords 1
for step $step do
  midplane outer
done
expressions $BOUNDARY
for step $step do
  boundary_quantities
done
EOF
  "$POSTPROC" < "tfit_pp_$step.pp" > "tfit_pp_$step.log" 2>&1
  local rc=$?
  local f1=postproc/exprs_outer-midplane_s$(printf %06d "$step").dat
  local f2=postproc/boundary_quantities_s$(printf %06d "$step").dat
  if [ $rc -ne 0 ] || [ ! -s "$f1" ] || [ ! -s "$f2" ]; then
    log "postproc failed for step $step (see tfit_pp_$step.log)"
    return 1
  fi
  return 0
}

crash() {                            # $1 = first step of the failed segment, $2 = reason
  local first=$1 n dir
  n=$(python3 -c "import json; print(json.load(open('$BASE/tfit_state.json'))['crashes'] + 1)")
  dir=tfit_crash_$n
  mkdir -p "$dir"
  for f in jorek[0-9]*.h5; do
    s=$(echo "$f" | sed -E 's/jorek0*([0-9]+)\.h5/\1/')
    [ "$s" -gt "$first" ] && mv "$f" "$dir/"
  done
  [ -f logfile ] && cp logfile "$dir/"
  cp "$NAMELIST" "$dir/"
  log "segment from step $first failed: $2. Step files written by it moved to $dir."
  log "If the job hit the wall time rather than crashing, raise WALLTIME or lower NSTEP_SEG and re-run."
  python3 "$PY" crash --cfg "$CFG" | tee -a "$LOG"
}

# ---------------------------------------------------------------- pre-flight
log "transport_fit.sh starting in $BASE (RUN_MODE=$RUN_MODE)"
for f in "$PY" "$NAMELIST"; do
  [ -f "$f" ] || { log "missing $f (copy util/transport_fit.py into $BASE)"; exit 1; }
done
if [ "$RUN_MODE" = toy ]; then
  [ -f "$TOY" ] || { log "missing $TOY"; exit 1; }
else
  [ -x "$POSTPROC" ] || { log "missing or not executable: $POSTPROC"; exit 1; }
  [ -x "$BASE/$EXE" ] || { log "missing or not executable: $BASE/$EXE"; exit 1; }
fi
write_cfg
eval "$(python3 "$PY" exprs --cfg "$CFG")"          # MIDPLANE=..., BOUNDARY=...
python3 "$PY" check --cfg "$CFG" | tee -a "$LOG"
[ "${PIPESTATUS[0]}" -eq 0 ] || { log "check failed, fix the problems above"; exit 1; }

if [ ! -f tfit_ref_midplane.dat ]; then
  log "exporting the reference step $REF_STEP"
  postproc_step "$REF_STEP" || exit 1
  python3 "$PY" reference --cfg "$CFG" --step "$REF_STEP" | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -eq 0 ] || { log "reference failed"; exit 1; }
fi

# ---------------------------------------------------------------- loop
while true; do
  eval "$(python3 "$PY" status --cfg "$CFG")"       # ITER NSTEP DONE CLAMP CRASHES RESTART_PARTICLES NOUT
  [ "$DONE" = 1 ] && break
  first=$(newest_step)
  [ -n "$first" ] || { log "no jorek??????.h5 step file to restart from"; exit 1; }
  cp "jorek$(printf %06d "$first").h5" jorek_restart.h5
  python3 "$PY" set-namelist --cfg "$CFG" --set "restart=.t." --set "nstep_n=$NSTEP" \
    --set "restart_particles=$RESTART_PARTICLES" || exit 1
  log "iteration $((ITER + 1)): segment of $NSTEP steps from step $first (clamp $CLAMP, crashes $CRASHES)"
  if ! run_segment; then
    crash "$first" "the run exited with an error"
    continue
  fi
  last=$(newest_step)
  if [ "$last" -le "$first" ]; then
    crash "$first" "no new step file was written"
    continue
  fi
  [ -f logfile ] && cp logfile "tfit_logfile_iter$(printf %02d $((ITER + 1)))"
  mid=$(( first + ( (last - first) / 2 / NOUT ) * NOUT ))
  [ "$mid" -gt "$first" ] || mid=$last
  postproc_step "$mid" || exit 1
  postproc_step "$last" || exit 1
  python3 "$PY" update --cfg "$CFG" --first "$first" --mid "$mid" --last "$last" | tee -a "$LOG"
  code=${PIPESTATUS[0]}
  case $code in
    0)  ;;
    10) log "converged"; break ;;
    20) log "state not steady: the segment is rerun with more steps" ;;
    30) log "stopped (iteration limit or repeated failures)"; break ;;
    *)  log "update failed with code $code"; exit 1 ;;
  esac
done

python3 "$PY" report --cfg "$CFG" | tee -a "$LOG"
log "done. Profile files are in place; originals kept as <file>.user, per-iteration backups as <file>.iterNN"
