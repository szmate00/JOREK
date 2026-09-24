#!/bin/bash
# Extract the floating_u wall diagnostics for one boundary type from a JOREK log into a text file:
#   1. per-step table for the type: step, t, min rho (R,Z), vE.n in (R,Z), vE.n out, min Te (R,Z), inflow, exb>cs,
#      sink/Bohm, and from [mach1]: |res| max, max M (R,Z)
#   2. the first step where min rho of the type is negative (and where min Te is)
#   3. the [wall prof] lines of the type inside an R window, for the last N steps that have a profile
# usage: extract_wall_diag.sh LOGFILE [OUT=wall_diag_extract.txt] [TYPE=1] [RMIN=1.585] [RMAX=1.600] [NLAST=20]
set -e
log=${1:?logfile}; out=${2:-wall_diag_extract.txt}; typ=${3:-1}; rmin=${4:-1.585}; rmax=${5:-1.600}; nlast=${6:-20}

{
echo "# source: $log   type: $typ   profile window R in [$rmin, $rmax]   last $nlast profile steps"
echo
echo "## 1. per-step table, boundary type $typ"
printf "%8s %12s | %11s %8s %8s | %11s %8s %8s | %11s | %11s %8s %8s | %7s %7s %10s | %11s %8s %8s %8s\n" \
  step t_now min_rho R Z vEn_in R Z vEn_out min_Te R Z inflow "exb>cs" sink/Bohm res_max max_M R Z
awk -v T="$typ" '
  $1=="[floating_u]" && $2=="step" { step=$3; t=$5; next }
  $1=="[floating_u]" && $2==T && step!="" {
    fu[step]=sprintf("%11s %8s %8s | %11s %8s %8s | %11s | %11s %8s %8s | %7s %7s %10s", $12,$13,$14, $9,$10,$11, $6, $15,$16,$17, $18,$19,$20)
    tt[step]=t; if (!(step in seen)) { order[++n]=step; seen[step]=1 } ; next }
  $1=="[mach1]" && $2==T && step!="" { m1[step]=sprintf("%11s %8s %8s %8s", $5, $8, $9, $10); next }
  END { for (i=1;i<=n;i++) { s=order[i]; printf "%8s %12s | %s | %s\n", s, tt[s], fu[s], m1[s] } }
' "$log"
echo
echo "## 2. first steps with a non-positive minimum on type $typ"
awk -v T="$typ" '
  $1=="[floating_u]" && $2=="step" { step=$3; next }
  $1=="[floating_u]" && $2==T && step!="" {
    if (!r && $12+0 <= 0) { r=1; printf "first min rho <= 0 : step %s  rho %s at (%s, %s)\n", step, $12, $13, $14 }
    if (!e && $15+0 <= 0) { e=1; printf "first min Te  <= 0 : step %s  Te  %s at (%s, %s)\n", step, $15, $16, $17 }
    if (!f && $3+0  > 1e-3) { f=1; printf "first |u-uf| > 1 mV (floor branch flip or row failure): step %s  %s V at (%s, %s)\n", step, $3, $4, $5 } }
  END { if (!r) print "min rho never <= 0"; if (!e) print "min Te never <= 0"; if (!f) print "|u-uf| never > 1 mV" }
' "$log"
echo
echo "## 3. [wall prof] type $typ, R in [$rmin, $rmax], last $nlast profile steps"
echo "# step  type  R  Z  rho[1e20]  Ti[eV]  Te[eV]  Phi[V]  Vpar*Bn  vE.n  cs|b.n|  vn[m/s]  b.n  M  dTe/ds[eV/m]  dPhi/ds[V/m]"
awk -v T="$typ" -v A="$rmin" -v B="$rmax" -v N="$nlast" '
  $1=="[wall" && $2=="prof]" && $3=="step" { step=$4; if (!(step in seen)) { order[++n]=step; seen[step]=1 } ; next }
  $1=="[wall" && $2=="prof]" && $3==T && step!="" && $4+0>=A && $4+0<=B { line[step]=line[step] sprintf("%8s  %s\n", step, substr($0, index($0,$3))) ; next }
  END { i0 = (n>N) ? n-N+1 : 1
        for (i=i0;i<=n;i++) { s=order[i]; if (s in line) printf "%s", line[s] }
        if (n==0) print "# no [wall prof] output in this log (set wall_diag_profile_every > 0)" }
' "$log"
} > "$out"

echo "wrote $out ($(wc -l < "$out") lines)"
