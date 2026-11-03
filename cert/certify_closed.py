# Rigorous lower bound for C16 (16a) using the FULL infinite trial function
#   u = sum_{k in Z^2 \ 0} c_k e^{ik.x},  c_k = 1/(n(1+mu n)),  n=|k|^2,  mu = MUNUM/MUDEN (exact rational)
# on [0,2pi]^2. With S0=sum c_k, S1=sum n c_k^2, S2=sum n^2 c_k^2:
#   ||grad u||^2 = 4pi^2 S1, ||Lap u||^2 = 4pi^2 S2, u(0)=S0 <= ||u||_inf, delta = S2/S1,
#   L >= S0^2/(pi S1) - ln(delta) - ln(1+ln delta)   (RHS decreasing in delta, delta>=1).
# All infinite sums are evaluated with NO truncation of the inner (y) direction: with lam=1/mu,
#   1/(n(1+mu n)) = 1/n - 1/(n+lam)
#   n/(n^2(1+mu n)^2) = 1/n - 1/(n+lam) - lam/(n+lam)^2
#   n^2/(n^2(1+mu n)^2) = lam^2/(n+lam)^2
# and the exact row sums over y in Z (a>0):
#   A(a) = sum_y 1/(y^2+a^2)   = (pi/a) coth(pi a)
#   B(a) = sum_y 1/(y^2+a^2)^2 = pi/(2a^3) coth(pi a) + pi^2/(2a^2) csch(pi a)^2
# Row x=0 (y != 0) uses zeta(2)=pi^2/6 per side. Rows 1<=|x|<=X are summed term by term in arb
# (ball arithmetic, every operation outward-rounded). Rows |x|>X: f0(x)=pi*lam/(x b (x+b)) and
# h(x)=pi/(2 b^3), b=sqrt(x^2+lam), are positive decreasing, so sum_{x>X} f in [int_{X+1}^inf f, int_X^inf f]
# (closed forms below); the hyperbolic corrections coth-1, csch^2 for x>X are < 10*e^{-2 pi X}
# per unit and are covered by an explicit added radius.
import sys
from flint import arb, ctx
ctx.prec = 128
MUNUM = int(sys.argv[1]); MUDEN = int(sys.argv[2]); X = int(sys.argv[3])
mu = arb(MUNUM) / MUDEN
lam = arb(MUDEN) / MUNUM
pi = arb.pi()
def A(a):
    return pi / a * (pi * a).coth()
def B(a):
    s = (pi * a).csch()
    return pi / (2 * a**3) * (pi * a).coth() + pi**2 / (2 * a**2) * s * s
# x = 0 row, y != 0
c = lam.sqrt()
z2 = pi**2 / 3                      # sum_{y!=0} 1/y^2
t0 = A(c) - 1 / lam                 # sum_{y!=0} 1/(y^2+lam)
t2 = B(c) - 1 / lam**2              # sum_{y!=0} 1/(y^2+lam)^2
R0 = [z2 - t0, z2 - t0 - lam * t2, lam**2 * t2]
# rows 1..X
acc = [arb(0), arb(0), arb(0)]
for x in range(1, X + 1):
    xa = arb(x)
    b = (xa * xa + lam).sqrt()
    Ax = A(xa); Ab = A(b); Bb = B(b)
    d = Ax - Ab
    acc[0] += d
    acc[1] += d - lam * Bb
    acc[2] += lam**2 * Bb
# tails x > X (one side)
def F0(Y):   # int_Y^inf (1/x - 1/sqrt(x^2+lam)) dx = ln((Y+b)/(2Y))
    Y = arb(Y); b = (Y * Y + lam).sqrt(); return ((Y + b) / (2 * Y)).log()
def H(Y):    # int_Y^inf (x^2+lam)^{-3/2} dx = (1 - Y/b)/lam
    Y = arb(Y); b = (Y * Y + lam).sqrt(); return (1 - Y / b) / lam
def encl(lo, hi):
    return arb((lo + hi) / 2) + arb(0, 1) * ((hi - lo) / 2).upper() if False else (lo.union(hi))
T0 = (pi * F0(X + 1)).union(pi * F0(X))            # sum_{x>X} pi(1/x - 1/b)
TH = (pi / 2 * H(X + 1)).union(pi / 2 * H(X))      # sum_{x>X} pi/(2 b^3)
hyp = arb(0, arb(10) ** -200) if X >= 100 else None  # coth/csch corrections, < 1e3*e^{-2piX} total
assert hyp is not None
tails = [T0 + hyp, T0 - lam * TH + hyp, lam**2 * TH + hyp]
S = [R0[j] + 2 * (acc[j] + tails[j]) for j in range(3)]
S0, S1, S2 = S
delta = S2 / S1
assert delta > 1
val = S0 * S0 / (pi * S1) - delta.log() - (1 + delta.log()).log()
print("mu", f"{MUNUM}/{MUDEN}", "X", X)
for j in range(3): print(f"S{j}", S[j].str(25))
print("delta", delta.str(20))
print("L in", val.str(20))
print("L >=", val.lower().str(18))
