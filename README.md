<p align="center">
  <img alt="An open specimen cabinet" src=".github/assets/specimen-cabinet-dark.jpg" width="420">
</p>

<h1 align="center">Conjectures</h1>

<p align="center">
  Open conjectures tested by Lugel Labs — each with its paper, an exact verifier and the output to expect.
</p>

<p align="center">
  <a href="https://lugel.org/blog/math-conjectures"><img alt="Map" src="https://img.shields.io/badge/map-lugel.org-1f1f1f?style=flat-square"></a>
  <a href="#conjectures"><img alt="Conjectures" src="https://img.shields.io/badge/conjectures-1-7c6fd6?style=flat-square"></a>
  <img alt="Python" src="https://img.shields.io/badge/python-3.10%2B-3776ab?style=flat-square&logo=python&logoColor=white">
  <a href="LICENSE"><img alt="Code: MIT" src="https://img.shields.io/badge/code-MIT-2e7d32?style=flat-square"></a>
  <a href="LICENSE-CC-BY-4.0"><img alt="Papers: CC BY 4.0" src="https://img.shields.io/badge/papers-CC%20BY%204.0-ef9421?style=flat-square"></a>
</p>

---

<a id="conjectures"></a>

| Conjecture | Area | Date | Status | Links |
| :-- | :-- | :-- | :-: | :-: |
| **Binder–Vecchi real-rootedness** | Matroids | Nov 2026 | ![Disproved](https://img.shields.io/badge/disproved-b42318?style=flat-square) | [![Code](https://img.shields.io/badge/code-1f1f1f?style=flat-square&logo=github)](https://github.com/Lugel-Labs/conjectures/tree/binder-vecchi/binder-vecchi) [![Note](https://img.shields.io/badge/note-7c6fd6?style=flat-square)](https://lugel.org/blog/math-conjectures/binder-vecchi) |

## How this repository works

- **One branch per conjecture**, named after it. `main` holds only this index.
- **One folder per conjecture** on its branch, with the same name.
- **Every folder ships** the paper and its LaTeX source, a `requirements.txt`, a verifier, its expected output and a `SHA256SUMS.txt`.

```sh
git clone --branch binder-vecchi https://github.com/Lugel-Labs/conjectures.git
cd conjectures/binder-vecchi
```

> [!NOTE]
> The papers carry the proofs. The code checks specific computations; it is not a formal proof. Results are not independently reviewed unless a folder says so.

## License

Code is released under the [MIT License](LICENSE). Papers and their LaTeX sources are released under [CC BY 4.0](LICENSE-CC-BY-4.0).

---

<p align="center">
  <a href="https://r-0.org"><picture><source media="(prefers-color-scheme: dark)" srcset=".github/assets/ring0-wordmark-white.png"><img alt="Developed by Ring0 Labs" src=".github/assets/ring0-wordmark-black.png" height="36"></picture></a>
  &nbsp;&nbsp;&nbsp;
  <a href="https://lugel.org"><picture><source media="(prefers-color-scheme: dark)" srcset=".github/assets/lugel-lockup-white.png"><img alt="Published by Lugel Labs" src=".github/assets/lugel-lockup-black.png" height="36"></picture></a>
  <br>
  <sub>Developed by Ring0 Labs · Published by Lugel Labs</sub>
</p>
