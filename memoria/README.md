# X17 Typst report (`memoria/`)

This folder contains the final English Typst paper and the auxiliary material used by it.

## Files

- `paper.typ`: main Typst manuscript.
- `refs.bib`: bibliography used by the manuscript.
- `figures/geometry/`: geometry views copied from `geo/`.
- `figures/rootfits/`: ROOT fit PDFs copied from the available `results/<geometry>/plots_signal_background_fit/` folders.
- `figures/sampling/`: generator-level and acceptance PDFs produced by `scripts/make_memoria_figures.py`.
- `tables/`: CSV tables and Typst table fragments used by `paper.typ`.
- `scripts/make_memoria_figures.py`: pure-Python reproducible generator-level figures and independent geometrical acceptance cross-check.
- `scripts/extract_root_summaries.C`: ROOT macro to extract official `summary` and `fit_summary` trees into CSV/Typst tables when ROOT is available.

## Regenerate figures and tables

From the project root:

```bash
python3 memoria/scripts/make_memoria_figures.py
```

To extract the official ROOT summaries from the existing ROOT files:

```bash
root -l -q 'memoria/scripts/extract_root_summaries.C()'
```

This writes:

```text
memoria/tables/root_geometrical_efficiencies.csv
memoria/tables/root_geometrical_efficiencies.typ
memoria/tables/root_fit_summary.csv
memoria/tables/root_fit_summary.typ
```

The current manuscript uses the deterministic Python-generated table fragments by default because they can be regenerated without ROOT. To use ROOT-derived tables instead, replace the `input("tables/geometrical_efficiencies.typ")` and `input("tables/reconstructed_yields.typ")` lines in `paper.typ` with the corresponding ROOT table fragments.

## Compile

```bash
typst compile memoria/paper.typ memoria/paper.pdf
```

or from inside `memoria/`:

```bash
typst compile paper.typ paper.pdf
```

All generated figures are PDFs, matching the project convention.
