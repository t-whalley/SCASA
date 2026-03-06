"""
Regression tests for SCASA shape complementarity against TCR3d reference values.

Reference SC values are taken from the TCR3d database (test/data/values/tcr3d.csv),
which reports scores computed by the CCP4 SC program using MSMS-based surface
generation. SCASA uses a ConvexHull triangulation approach, which produces
systematically higher scores by approximately 0.05–0.15 due to the convex surface
approximation being unable to represent concave interface regions.

Chain assignments (TCR alpha/beta vs antigen/MHC) are taken from STCRDab
(test/data/values/stcrdab.tsv).

Each test:
  - Requires the PDB file to exist at test/data/<PDBID>.pdb (lowercase)
  - Skips automatically if the PDB file is not present
  - Skips if the PDB ID is not found in STCRDab (no chain information)
  - Passes if |SCASA_SC - reference_SC| <= TOLERANCE

The tolerance of 0.15 accounts for the expected systematic offset between
SCASA and CCP4 SC, plus run-to-run stochasticity from random surface sampling.
To tighten this, increase dot density via the density parameter.
"""

import csv
import os
import pytest

# Absolute tolerance between SCASA SC and CCP4 SC reference values.
# SCASA scores are typically 0.05–0.15 higher than CCP4 SC due to the
# ConvexHull surface approximation.
TOLERANCE = 0.15

DATA_DIR    = os.path.join(os.path.dirname(__file__), "data")
VALUES_DIR  = os.path.join(DATA_DIR, "values")
TCR3D_CSV   = os.path.join(VALUES_DIR, "tcr3d.csv")
STCRDAB_TSV = os.path.join(VALUES_DIR, "stcrdab.tsv")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_reference_values():
    """
    Load reference SC values from TCR3d CSV.
    Returns dict mapping lowercase PDB ID -> float SC value.
    Entries with missing SC values are excluded.
    """
    values = {}
    with open(TCR3D_CSV, newline="") as f:
        for row in csv.DictReader(f):
            pdb = row["PDB ID"].strip().lower()
            sc_str = row["Shape complementarity"].strip()
            if sc_str:
                try:
                    values[pdb] = float(sc_str)
                except ValueError:
                    pass
    return values


def load_chain_assignments():
    """
    Load chain assignments from STCRDab TSV.
    Returns dict mapping lowercase PDB ID -> dict with keys:
        tcr_chains  : str  e.g. "DE"   (alpha + beta, or gamma + delta)
        other_chains: str  e.g. "ABC"  (antigen + MHC chain 1 + MHC chain 2)
    Only entries where both TCR chains AND at least one opposing chain are
    present are included. NA values are ignored.
    """
    assignments = {}
    with open(STCRDAB_TSV, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pdb = row["pdb"].strip().lower()
            tcr_type = row["TCRtype"].strip()

            # Resolve TCR chains based on receptor type
            if tcr_type == "abTCR":
                tcr_chains = [row["Achain"], row["Bchain"]]
            elif tcr_type == "gdTCR":
                tcr_chains = [row["Dchain"], row["Gchain"]]
            else:
                continue  # skip NA / unknown TCR types

            # Filter out NA values
            tcr_chains = [c for c in tcr_chains if c and c != "NA"]
            if len(tcr_chains) < 2:
                continue

            # Collect opposing chains: antigen + MHC subunits, filtering NAs
            other_chains = []
            for key in ("antigen_chain", "mhc_chain1", "mhc_chain2"):
                val = row[key].strip()
                if val and val != "NA":
                    other_chains.append(val)

            if not other_chains:
                continue

            assignments[pdb] = {
                "tcr_chains":   "".join(tcr_chains),
                "other_chains": "".join(other_chains),
            }

    return assignments


# ---------------------------------------------------------------------------
# Build test parameters at collection time
# ---------------------------------------------------------------------------

_reference  = load_reference_values()
_chains     = load_chain_assignments()

def _make_params():
    params = []
    for pdb_id, ref_sc in _reference.items():
        # PDB files may be stored with uppercase IDs (e.g. 1FYT.pdb)
        pdb_file_upper = os.path.join(DATA_DIR, f"{pdb_id.upper()}.pdb")
        pdb_file_lower = os.path.join(DATA_DIR, f"{pdb_id.lower()}.pdb")
        pdb_file = pdb_file_upper if os.path.isfile(pdb_file_upper) else pdb_file_lower
        chain_info = _chains.get(pdb_id)
        params.append(
            pytest.param(
                pdb_id,
                pdb_file,
                chain_info,
                ref_sc,
                id=pdb_id.upper(),
            )
        )
    return params


# ---------------------------------------------------------------------------
# Regression test
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("pdb_id, pdb_file, chain_info, ref_sc", _make_params())
def test_sc_regression(pdb_id, pdb_file, chain_info, ref_sc):
    """
    Calculate SC for a TCR–pMHC complex and compare to the TCR3d reference value.
    """
    if not os.path.isfile(pdb_file):
        pytest.skip(f"PDB file not present: {pdb_file}")

    if chain_info is None:
        pytest.skip(f"{pdb_id.upper()} not found in STCRDab — no chain assignment available")

    # Import here so collection does not depend on the package being installed
    from scasa.scasa import Complex

    tcr    = chain_info["tcr_chains"]
    other  = chain_info["other_chains"]

    complex_obj = Complex(
        pdb_file=pdb_file,
        complex_1=tcr,
        complex_2=other,
        verbose=False,
    )
    sc = complex_obj.sc()

    assert sc is not None, f"sc() returned None for {pdb_id.upper()}"
    assert abs(sc - ref_sc) <= TOLERANCE, (
        f"{pdb_id.upper()}: SCASA SC={sc:.4f}, reference={ref_sc:.4f}, "
        f"diff={abs(sc - ref_sc):.4f} > tolerance={TOLERANCE} "
        f"(TCR={tcr}, other={other})"
    )
