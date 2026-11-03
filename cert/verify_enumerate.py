# Second, separate check of the 16a bound: direct lattice enumeration in sectors plus unit-square tail bounds in arb.
# (sector enumeration + unit-square tail bounds). Changes: large N, and a per-term margin derived
# from an explicit rounding count instead of 3*2^-40.
# T^2=[0,2pi]^2, u = sum_{k!=0} c_k e^{ik.x}, c_k = 1/(n(1+mu n)), n=|k|^2.
# ||grad u||^2 = 4pi^2 S1, ||Lap u||^2 = 4pi^2 S2, u(0)=S0;  4pi|u(0)|^2/||grad u||^2 = S0^2/(pi S1).
# L >= S0^2/(pi S1) - ln d - ln(1+ln d), d=S2/S1; RHS decreasing in d -> use upper bound of d.
# Rounding count (u=2^-53): n exact (<2^53); mu*n exact (mu dyadic, numerator*n < 2^53 asserted);
# d=1+mu n: 1 rounding; n*d: 2; g0=1/(n d): 3; g1=g0*g0*n: 8; g2=1/(d*d): 4 roundings (relative,
# first order). Terms positive, so a sum of terms each within rel e is within rel e.
# Row fsum: +1u, outer fsum: +1u, *4 exact. Total < 11u (+O(u^2)). Margin used: 2^-48 = 32u.
import sys, math, time
import numpy as np
from flint import arb, ctx
ctx.prec = 200
N = int(sys.argv[1]); MUNUM = int(sys.argv[2]); MUDEN = int(sys.argv[3])
assert MUDEN & (MUDEN - 1) == 0 and MUNUM.bit_length() + (2 * N * N).bit_length() <= 53
mu_f = MUNUM / MUDEN
mu = arb(MUNUM) / MUDEN
M = N * N
t0 = time.time()
rows = [[], [], []]; npts = 0
for x in range(1, N + 1):
    ymax = math.isqrt(M - x * x)
    y = np.arange(0, ymax + 1, dtype=np.float64)
    n = x * x + y * y
    d = 1.0 + mu_f * n
    g0 = 1.0 / (n * d)
    g1 = g0 * g0 * n
    g2 = 1.0 / (d * d)
    rows[0].append(math.fsum(g0)); rows[1].append(math.fsum(g1)); rows[2].append(math.fsum(g2))
    npts += len(y)
eps = arb(2) ** -48
fin = [arb(4 * math.fsum(rows[j])) * (1 + arb(0, 1) * eps) for j in range(3)]
a = arb(2).sqrt() / 2; pi = arb.pi()
def Isg(j, B):
    U = B * B
    if j == 0: return (1 + 1 / (mu * U)).log() / 2
    if j == 1: return ((1 + 1 / (mu * U)).log() - 1 / (1 + mu * U)) / 2
    return 1 / (2 * mu * (1 + mu * U))
def Ig_up(j, B): return Isg(j, B) / B
lo = []; hi = []
for j in range(3):
    Bu = N - 2 * a; Bl = N + 2 * a
    tu = 2 * pi * (Isg(j, Bu) + a * Ig_up(j, Bu))
    tl = 2 * pi * (Isg(j, Bl) - a * Ig_up(j, Bl))
    lo.append(arb((fin[j] + tl).lower())); hi.append(arb((fin[j] + tu).upper()))
d_up = arb((hi[2] / lo[1]).upper())
assert d_up > 1
val = lo[0] * lo[0] / (pi * hi[1]) - d_up.log() - (1 + d_up.log()).log()
d_lo = arb((lo[2] / hi[1]).lower())
val_up = hi[0] * hi[0] / (pi * lo[1]) - d_lo.log() - (1 + d_lo.log()).log()
print("N", N, "mu", f"{MUNUM}/{MUDEN}", "points", 4 * npts + 0, "secs", round(time.time() - t0, 1))
for j in range(3): print(f"S{j} in [{lo[j].mid().str(20)}, {hi[j].mid().str(20)}]")
print("delta <=", d_up.str(17))
print("L >= ", val.lower().str(18))
print("(value of this u is <= ", val_up.upper().str(18), ")")
