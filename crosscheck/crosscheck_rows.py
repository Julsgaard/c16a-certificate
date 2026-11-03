# Cross-check the closed-form row sums A(a), B(a) against direct summation (|y|<=Y plus 1/Y tail), and
# the full S_j at mu=128031/2^20 against an earlier, separately computed N=4000 enclosure.
import math
from flint import arb
pi=arb.pi()
for a in [0.7,1.0,2.5,9.3]:
    aa=arb(a)
    A=pi/aa*(pi*aa).coth(); s=(pi*aa).csch(); B=pi/(2*aa**3)*(pi*aa).coth()+pi**2/(2*aa**2)*s*s
    Y=10**6
    dA=1/a**2+2*math.fsum(1/(y*y+a*a) for y in range(1,Y+1))+2/Y
    dB=1/a**4+2*math.fsum(1/(y*y+a*a)**2 for y in range(1,Y+1))
    print(a,"A",A.mid().str(16),dA,"B",B.mid().str(16),dB)
