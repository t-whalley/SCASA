# SCASA: Shape Complementarity and Available Surface Area

SCASA is a Python package for calculating two geometric properties of protein complexes from PDB files: Shape Complementarity (SC) and Available/Buried Surface Area (ASA/BSA). It can be used as a standalone command-line tool or imported as a Python module.

---

## What is Shape Complementarity?

Shape Complementarity (SC) [(Lawrence and Colman, 1993)](https://pubmed.ncbi.nlm.nih.gov/8263940/) is a measure of the geometric "goodness of fit" between two protein interfaces. It ranges from 0 to 1, where values closer to 1 indicate tightly complementary surfaces. It has been widely used to characterise the quality of antibody/antigen and T-cell receptor/antigen interfaces.

SC is computed by:
1. Selecting the interface atoms from each complex (those within a distance threshold of the opposing surface)
2. Generating surface dots by randomly sampling a ConvexHull triangulation of each interface
3. Estimating a surface normal at each dot via PCA on its 10 nearest dot-cloud neighbours
4. For each dot on surface A, finding the nearest dot on surface B and computing the dot product of their normals
5. The final SC score is the mean of the medians of S(A→B) and S(B→A)

### Comparison with CCP4 SC

SCASA implements the same Connolly molecular surface algorithm as CCP4 SC, translated directly from the original Fortran `mds` subroutine (Copyright Michael Connolly, 1986). The surface generation produces the same three types of surface dots:

- **Convex** — contact surface on each atom's VdW shell
- **Toroidal** — probe rolling between two adjacent atoms
- **Concave** — re-entrant patch where the probe nestles between three atoms

SCASA uses the same radii file (`sc_radii.lib`), probe radius (1.7 Å), dot density (15/Å²), and weight factor (0.5) as CCP4 SC. In place of the CCP4 trim band (which removes peripheral buried dots near accessible dots), SCASA uses an equivalent 1.5 Å inter-surface distance filter, which produces the same practical effect.

At default settings, SCASA gives **SC ≈ 0.54** for 1FYT vs **0.56** from CCP4 SC — a difference of ~0.02, within the expected variation from floating-point differences and random surface sampling. Scores are directly comparable to CCP4 SC values.

To match CCP4 SC defaults exactly, use:

```bash
SCASA sc --pdb structure.pdb --complex_1 DE --complex_2 ABC --dot-density 15
```

## What is (Buried or Available) Surface Area?

Buried Surface Area (BSA; reviewed by [Ali et al., 2014](https://pubmed.ncbi.nlm.nih.gov/24678666/)) measures the total surface area (in Å²) of a complex that becomes buried upon binding to another complex. Available Surface Area (ASA) is the reciprocal — the total surface area remaining accessible to solvent.

Accessibility is defined relative to the solvent-accessible surface area (SASA), which imagines a probe sphere (representing a solvent molecule) rolling along the protein surface. Any region where the probe can pass unimpeded is considered accessible.

SCASA calculates ASA using the Shrake-Rupley algorithm implemented in Biopython. BSA is then derived by comparing the ASA of each complex in isolation against its ASA when bound.

---

## Installation

Requires Python ≥ 3.9.

```bash
pip install .
```

---

## Usage

### Command-line tool

SCASA provides two subcommands: `sc` for shape complementarity and `asa` for surface area.

#### Common arguments (both subcommands)

| Flag | Short | Required | Description |
|------|-------|----------|-------------|
| `--pdb` | `-P` | Yes | Path to PDB file of the complex |
| `--complex_1` | `-C1` | Yes | Chains of the first complex (e.g. `DE` or `ABC`). Multiple chains must be supplied as a single concatenated string |
| `--complex_2` | `-C2` | No | Chains of the second complex. Defaults to all remaining chains in the PDB file |
| `--verbose` | `-v` | No | Print additional progress messages |

#### `SCASA sc` — Shape Complementarity

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--distance` | `-D` | `8.0` | Interface cutoff in Å. Atoms with no neighbour within this distance of the opposing surface are excluded |
| `--dot-density` | `-Dd` | `1.5` | Surface dot sampling density (dots per Å² of interface area) |
| `--plot` | `-pl` | — | Generate a histogram plot of the SC function distribution |

Example:
```bash
SCASA sc --pdb test/data/1FYT.pdb --complex_1 DE --complex_2 ABC
```

#### `SCASA asa` — Available/Buried Surface Area

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--level` | `-L` | `R` | Granularity of output: `S` (whole complex), `C` (per chain), `R` (per residue), or `A` (per atom) |

Example:
```bash
SCASA asa --pdb test/data/1FYT.pdb --complex_1 DE --complex_2 ABC --level R
```

---

### Python module

SCASA can be imported and used directly in Python via the `Complex` class:

```python
from scasa.scasa import Complex

complex = Complex(
    pdb_file="test/data/1FYT.pdb",
    complex_1="DE",
    complex_2="ABC",   # optional — defaults to all remaining chains
    distance=8.0,      # interface cutoff in Å
    density=1.5,       # dot sampling density per Å²
    verbose=True,
)

# Calculate shape complementarity
complex.sc()

# Calculate ASA/BSA (requires sub-PDB files written to tmp/)
complex.create_sub_pdbs()
complex.complex_sasa()
```

---

## Interpreting SC scores

As a rough guide based on published literature:

| SC score (SCASA) | Interpretation |
|------------------|---------------|
| 0.70 – 0.80 | Tightly complementary (e.g. antibody–antigen) |
| 0.60 – 0.70 | Typical protein–protein interface |
| 0.45 – 0.60 | Loosely packed or transient complex |
| < 0.45 | Poor fit; may indicate a crystal contact rather than a biological interface |

Note: these ranges apply when using the default density of 15 dots/Å² (--dot-density 15), which matches CCP4 SC. At the lower default density of 1.5, scores will be somewhat lower.

---

## Contact

[Tom Whalley](mailto:whalleyt@cardiff.ac.uk)
