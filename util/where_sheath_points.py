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

step, a = blocks[-1]
d = {n: a[:, i] for i, n in enumerate(names)}
R, Z = d['R'], d['Z']
bn = np.abs(d['Bnorm']) / d['B_abs']
ratio = np.abs(d['Jpar']) / np.abs(d['Jpar_ionsat'])
ag_all = np.sign(d['Jpar']) == np.sign(d['Jpar_ionsat'])

print(f"--- step {step} ---")
k = np.ones(len(R), bool)
for lab, c in [('ALL boundary pts', k.copy()), ('|B.n|/B > 0.02', bn > 0.02),
               ('ne > 1e18', d['ne'] > 1e18), ('T_e > 2 eV', d['T_e'] > 2.0),
               ('ratio 0.3..3', (ratio > 0.3) & (ratio < 3))]:
    k &= c
    if k.sum():
        print(f"  {lab:<18s} {k.sum():5d}   R {R[k].min():.3f}..{R[k].max():.3f}"
              f"   Z {Z[k].min():+.3f}..{Z[k].max():+.3f}")
    else:
        print(f"  {lab:<18s} {k.sum():5d}")

print("\n  R distribution of the WHOLE wall (is there anything out at large R at all?)")
h, e = np.histogram(R, bins=10)
for i in range(10):
    print(f"    {e[i]:.3f}-{e[i+1]:.3f} : {h[i]:5d} {'#'*int(round(40*h[i]/max(h)))}")

Rk, agk = R[k], ag_all[k]
print("\n  agreement vs R over the SURVIVORS (two clusters = two targets)")
e = np.linspace(Rk.min(), Rk.max(), 9)
for i in range(8):
    m = (Rk >= e[i]) & ((Rk <= e[i+1]) if i == 7 else (Rk < e[i+1]))
    if m.sum():
        print(f"    {e[i]:.3f}-{e[i+1]:.3f} : {m.sum():4d} pts, agree {100*agk[m].mean():5.1f} %"
              f"   Z {Z[k][m].min():+.3f}..{Z[k][m].max():+.3f}")
    else:
        print(f"    {e[i]:.3f}-{e[i+1]:.3f} :    0 pts")
