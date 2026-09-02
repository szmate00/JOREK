import re, sys, numpy as np
BLOCK = re.compile(r"time step #\s*(\d+).*?t_now\s*=\s*([-\dEe.+]+)")
def isnum(t):
    try: float(t); return True
    except ValueError: return False
names, rows, step, blocks = None, [], None, []
for line in open(sys.argv[1]):
    s = line.strip()
    if not s: continue
    if s.startswith('#'):
        m = BLOCK.search(s)
        if m:
            if step is not None and rows: blocks.append((step, np.array(rows, float)))
            rows, step = [], int(m.group(1))
        else:
            c = s.lstrip('#').split()
            if c and not isnum(c[0]): names = c
        continue
    rows.append([float(v) for v in s.split()])
if step is not None and rows: blocks.append((step, np.array(rows, float)))

ZX = float(sys.argv[2]) if len(sys.argv) > 2 else -0.93   # X-point Z
step, a = blocks[-1]
d = {n: a[:, i] for i, n in enumerate(names)}
R, Z, ne = d['R'], d['Z'], d['ne']
bn = np.abs(d['Bnorm']) / d['B_abs']
ratio = np.abs(d['Jpar']) / np.abs(d['Jpar_ionsat'])
ag = np.sign(d['Jpar']) == np.sign(d['Jpar_ionsat'])
base = (bn > 0.02) & (d['T_e'] > 2.0) & (ratio > 0.3) & (ratio < 3) & (Z < ZX)

print(f"--- step {step};  DIVERTOR ONLY (Z < {ZX}), scanning the ne threshold ---")
print(f"{'ne_min':>8s} {'pts':>5s}  {'R range':>15s}   {'HFS R<1.42':>22s}   {'LFS R>1.42':>22s}")
for nemin in (1e18, 3e17, 1e17, 3e16, 1e16, 0.0):
    k = base & (ne > nemin)
    if not k.sum():
        print(f"{nemin:8.0e} {0:5d}"); continue
    Rk, agk = R[k], ag[k]
    out = f"{nemin:8.0e} {k.sum():5d}  {Rk.min():6.3f}..{Rk.max():6.3f}"
    for lab, m in (("hfs", Rk < 1.42), ("lfs", Rk >= 1.42)):
        out += (f"   {agk[m].sum():4d}/{m.sum():<4d} agree {100*agk[m].mean():5.1f}%"
                if m.sum() else f"   {'-- none --':>22s}")
    print(out)

k = base & (ne > 1e16)
if k.sum():
    print("\n  LFS (R > 1.42) divertor points at ne > 1e16, resolved:")
    Rk, Zk, agk, nek = R[k], Z[k], ag[k], ne[k]
    m = Rk >= 1.42
    if m.sum():
        print(f"    {m.sum()} pts  R {Rk[m].min():.3f}..{Rk[m].max():.3f}"
              f"  Z {Zk[m].min():+.3f}..{Zk[m].max():+.3f}"
              f"  ne {nek[m].min():.2e}..{nek[m].max():.2e}"
              f"  agree {100*agk[m].mean():.1f} %")
    else:
        print("    none - the LFS target has no points passing even ne > 1e16")
