# C_16 (Brezis–Gallouet–Wainger remainder, 2D torus): lower bound 2.156225528

Trial function on T^2 = [0, 2π]^2, with no lattice cut-off:

    u = sum_{k in Z^2, k != 0} e^{ik·x} / (n (1 + μ n)),   n = |k|^2,   μ = 128042 / 2^20

With S0 = sum c_k, S1 = sum n c_k^2, S2 = sum n^2 c_k^2 (c_k the coefficients above), u(0) = S0 ≤ ||u||_∞,
4π u(0)^2 / ||∇u||^2 = S0^2 / (π S1) and δ = S2 / S1 ≈ 3.92884, so

    L ≥ S0^2 / (π S1) − ln δ − ln(1 + ln δ)

the same bound formula as the existing lower-bound rows on the 16a page. The right side falls as δ rises,
so an upper enclosure of δ is used.

## Check

```
python cert/certify_closed.py 128042 1048576 1000000
python cert/verify_enumerate.py 80000 128042 1048576
```

Both need python-flint (arb ball arithmetic); the second also numpy.

- `certify_closed.py` (about 5 s) splits each term by partial fractions in n and sums each row in y in
  closed form: sum_y 1/(y^2 + a^2) = (π/a) coth(πa), and the matching formula for the square. Rows
  |x| ≤ 10^6 are summed in arb; rows beyond are enclosed by integral comparison. Output in
  `cert/cert_mu128042_X1e6.txt`: L in [2.1562255281542061 ± 6.1e-17].
- `verify_enumerate.py` (about 6 minutes) is a separate check by a different route: direct enumeration of
  the lattice to radius N in sectors, with each float term given an explicit error margin, plus unit-square
  tail bounds in arb. Output in `cert/verify_enumerate_N80000.txt`: L ≥ 2.15622552815389.
- `crosscheck/` compares the closed-form row sums with direct summation.

The printed bound 2.156225528 is the lower of the two, rounded down.

## How it was found

This is the conditional extremal family of [BDZ2013], as used by [Pat2026] and PR #270. The gain comes from
dropping the |k| ≤ 4000 cut-off and re-tuning μ. Neighbouring μ (128038, 128046) are lower by about 1e-11.

Prepared with Claude (Anthropic, Claude Opus 5.5).
