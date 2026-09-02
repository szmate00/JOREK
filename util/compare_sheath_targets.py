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
ZX = float(sys.argv[2]) if len(sys.argv) > 2 else -0.93

step, a = blocks[-1]
d = {n: a[:, i] for i, n in enumerate(names)}
R, Z = d['R'], d['Z']
bn = np.abs(d['Bnorm']) / d['B_abs']
base = (bn > 0.02) & (d['T_e'] > 2.0) & (Z < ZX) & (d['ne'] > 1e16)

print(f"--- step {step}: the two divertor targets, side by side ---")
print("KEY: 'flow sign' = sign(Jpar_ionsat * Bnorm). SAME sign at both targets means")
print("     both legs flow the same way relative to the wall (normal). OPPOSITE means")
print("     the HFS leg has REVERSED parallel flow - plasma leaving the target.\n")
for lab, m in (("HFS  R<1.42", base & (R < 1.42)), ("LFS  R>1.42", base & (R >= 1.42))):
    if not m.sum():
        print(f"  {lab}: no points"); continue
    jp, js, bnorm = d['Jpar'][m], d['Jpar_ionsat'][m], d['Bnorm'][m]
    flow = np.sign(js * bnorm)
    agree = np.sign(jp) == np.sign(js)
    print(f"  {lab}:  {m.sum():4d} pts   R {R[m].min():.3f}..{R[m].max():.3f}"
          f"   Z {Z[m].min():+.3f}..{Z[m].max():+.3f}")
    print(f"      ne      median {np.median(d['ne'][m]):.3e} m^-3")
    print(f"      T_e     median {np.median(d['T_e'][m]):8.3f} eV"
          f"      T_i median {np.median(d['T_i'][m]):8.3f} eV" if 'T_i' in d else "")
    print(f"      Phi     median {np.median(d['Phi'][m]):+8.3f} V"
          f"   (Lambda*Te = {3*np.median(d['T_e'][m]):+.2f} V)")
    print(f"      Jpar    median {np.median(jp):+.4e} A/m^2")
    print(f"      Jsat    median {np.median(js):+.4e} A/m^2")
    print(f"      sign(Jpar)==sign(Jsat): {100*agree.mean():5.1f} %")
    print(f"      flow sign: {100*(flow > 0).mean():5.1f} % positive,"
          f" {100*(flow < 0).mean():5.1f} % negative")
    print(f"      -> net current is {'ION-side (into wall)' if np.median(jp*bnorm) * np.median(js*bnorm) > 0 else 'ELECTRON-side (opposite the ion flow)'}\n")
