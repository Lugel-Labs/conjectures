# Explicit counterexample in degree 1

**Claim.** Let $M=\mathrm{Tr}^7(\mathcal M_{47}\oplus\mathcal M_{47}\oplus U_{1,1})$, the simple rank-7 matroid of Bérczi–Coron–Eur–Li (Theorem C), on 9,761,573 elements. For every

$$p \ge p_0 = 2{,}144{,}540{,}694{,}648{,}271{,}312,$$

adding $p$ parallel elements gives a loopless matroid $M^{(p)}$ whose refined polynomial $\underline{\mathrm H}^1_{M^{(p)}}$ has a non-real zero.

## Run

```sh
python certify_i1.py          # about 2 seconds; needs mpmath (installed with SymPy)
```

Compare with `expected_output.txt`.

## How it works

1. **Exact Chow polynomial** (`chow_bcel.py`). $\underline{\mathrm H}_M$ is computed from the flag $f$-vector of $\mathcal M_q$ (BCEL, Corollary 4.11) through the direct sum and the truncation, then the Feichtner–Yuzvinsky basis. `threshold_check.txt` shows it reproduces both BCEL thresholds: real-rooted at $q=43$, not at $q=47$; log-concave at $q=89$, not at $q=97$.
2. **Binder's formula** (Prop. 3.6) in degree 1: $\underline{\mathrm H}^1_{M^{(p)}} = p\,\underline{\mathrm H}_M + \underline{\mathrm H}^1_M$.
3. **Circle.** Centre $\approx -384.12 + 124.71i$, radius $31.1767$, above the real axis. The circle is sampled at 8,192 exact rational points with a Lipschitz bound in between. This gives $|\underline{\mathrm H}_M| \ge 3.3\times10^{36}$ on the circle and an exact winding number of 1.
4. **Bound on $\underline{\mathrm H}^1_M$, never computed.** The Bergman complex is shellable, so reducing its face ring by a linear system of parameters gives an Artinian algebra of dimension $F$ (the number of maximal chains of proper flats). $\mathrm{Tor}_1$ is a subquotient of $F^{\,n-r}$ copies, so the coefficients of $\underline{\mathrm H}^1_M$ are nonnegative, lie in degrees $0..6$ and sum to at most $(n-r)F$.
5. **Rouché.** For $p \ge p_0$, $|\underline{\mathrm H}^1_M| < p\,|\underline{\mathrm H}_M|$ on the circle, so $\underline{\mathrm H}^1_{M^{(p)}}$ has a zero inside it, hence a non-real zero.

> [!NOTE]
> All steps are exact integer or rational arithmetic. Floating point is used only to choose the circle. $p_0$ comes from a crude bound and is far from sharp. The certificate relies on the cited theorems (BCEL Corollary 4.11, Binder Prop. 3.6, Franz's isomorphism); it is not a formal proof.
