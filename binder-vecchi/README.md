<p align="center">
  <img alt="An open specimen cabinet" src="../.github/assets/specimen-cabinet-dark.jpg" width="420">
</p>

<h1 align="center">Binder–Vecchi real-rootedness</h1>

<p align="center">
  <em>A remark on the real-rootedness of refined Hodge–Poincaré polynomials of matroids</em><br>
  Luis Ángel Gonzalez · Ring0 Labs · AI Swarms
</p>

<p align="center">
  <img alt="Disproved" src="https://img.shields.io/badge/status-disproved-b42318?style=flat-square">
  <a href="paper.pdf"><img alt="Paper" src="https://img.shields.io/badge/paper-PDF-1f1f1f?style=flat-square"></a>
  <a href="https://lugel.org/blog/math-conjectures/binder-vecchi"><img alt="Note" src="https://img.shields.io/badge/note-lugel.org-7c6fd6?style=flat-square"></a>
  <img alt="Python" src="https://img.shields.io/badge/python-3.10%2B-3776ab?style=flat-square&logo=python&logoColor=white">
  <a href="../LICENSE"><img alt="Code: MIT" src="https://img.shields.io/badge/code-MIT-2e7d32?style=flat-square"></a>
  <a href="../LICENSE-CC-BY-4.0"><img alt="Paper: CC BY 4.0" src="https://img.shields.io/badge/paper-CC%20BY%204.0-ef9421?style=flat-square"></a>
</p>

---

## The result

Binder and Vecchi conjectured that every refined Hodge–Poincaré polynomial $\underline{\mathrm H}^i_M(x)$ of a loopless matroid has only real zeros ([arXiv:2609.15946](https://arxiv.org/abs/2609.15946), Conj. 1.5).

**It fails in every homological degree.** For each $i \ge 0$ some loopless matroid has a polynomial in degree $i$ with a non-real zero.

- **$i = 0$** is already known: it is the Chow polynomial, disproved by Bérczi–Coron–Eur–Li ([arXiv:2610.17135](https://arxiv.org/abs/2610.17135)). Not claimed here.
- **$i \ge 1$** is new: adding parallel elements (Binder's formula, [arXiv:2412.05732](https://arxiv.org/abs/2412.05732v2)) and Rouché's theorem carry the failure to every degree.

## Explicit counterexample in degree 1

For $M=\mathrm{Tr}^7(\mathcal M_{47}\oplus\mathcal M_{47}\oplus U_{1,1})$ with **$p \ge 2{,}144{,}540{,}694{,}648{,}271{,}312$** parallel elements added, $\underline{\mathrm H}^1$ has a non-real zero. The exact certificate is in [`certificate/`](certificate) and runs in seconds.

## Reproduce

Python 3.10+, NumPy and SymPy. From this folder:

```sh
python -m pip install -r requirements.txt
python verify_results.py binder_vecchi_verificacion.zip
```

> [!IMPORTANT]
> Keep `binder_vecchi_verificacion.zip` unchanged: the verifier checks its SHA-256 before importing it. The PG(3,2) stage needs substantial time and memory.

### Expected output

| Matroid | Discriminant degree in $t$ | Smallest coefficient |
| :-- | --: | --: |
| AG(3,2) | 14 | 2,816,000 |
| M(W₄) | 16 | 70,848 |
| M(K₅) | 22 | 1,562,500 |

For PG(3,2), coefficients in increasing powers of $x$:

| Column | Coefficients |
| :-- | :-- |
| $\underline{\mathrm H}^9$ | `[0, 35, 2380, 1759]` |
| $\underline{\mathrm H}^{10}$ | `[0, 0, 245, 516]` |
| $\underline{\mathrm H}^{11}$ | `[0, 0, 0, 64]` |

Adjacent boundary ranks at exterior degree 4: 1,001 and 45,439. Discriminant of $\underline{\mathrm H}^9/x$: 5,418,140.

<details>
<summary><strong>Supplementary and exploratory scripts</strong></summary>

The verifier above is enough. These give alternative elimination or diagnostics:

```sh
python supplementary/check_cells.py binder_vecchi_verificacion.zip
python supplementary/pg_rational.py binder_vecchi_verificacion.zip
python supplementary/panel_rank.py binder_vecchi_verificacion.zip
```

Exploratory experiments need the archive extracted first:

```sh
python -m zipfile -e binder_vecchi_verificacion.zip experiments
python experiments/run_all.py          # add --full for the randomized modular run
python experiments/check_D.py
python experiments/check_perm.py
python experiments/scan_uniform.py
```

Modular ranks and numerical root scans are not rational certificates; the main verifier uses exact elimination over $\mathbb Q$. Sampled parameters do not prove a universal statement.

</details>

## Files

| File | What it is |
| :-- | :-- |
| `paper.pdf` · `paper.tex` · `references.bib` | The manuscript and its LaTeX source |
| `verify_results.py` | Exact verifier over $\mathbb Q$ |
| `binder_vecchi_verificacion.zip` | Sealed computation modules, checked by SHA-256 |
| `certificate/` | Explicit degree-1 counterexample with an exact certificate |
| `verification_smoke_test.txt` | Recorded output of a release run |
| `supplementary/` · `experiments/` | Alternative and exploratory checks |
| `SHA256SUMS.txt` | Checksums of every file here |
| `CITATION.cff` | How to cite this work |

## Scope

> [!NOTE]
> - The proof of Theorem 3 is in the paper. The software is **not** a formal proof, nor a reconstruction of the rank-seven seed counterexample.
> - The release run reproduced the three rank-four discriminants (Python 3.13.5, NumPy 2.3.5, SymPy 1.14.0). The PG(3,2) stage was not completed in that time-limited run; see `verification_smoke_test.txt`.
> - Supplementary and exploratory scripts were not all rerun for this package.
> - Priority and independent review have not been established by these checks.

## License

Code: [MIT](../LICENSE). Paper and LaTeX source: [CC BY 4.0](../LICENSE-CC-BY-4.0).

---

<p align="center">
  <sub>Developed by <a href="https://r-0.org">Ring0 Labs</a> · Published by <a href="https://lugel.org">Lugel Labs</a> · <a href="../README.md">All conjectures</a></sub>
</p>
