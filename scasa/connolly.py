"""
Connolly molecular dot surface (MDS) — Python translation of sc.f mds subroutine.

Translated from Fortran by Michael Connolly (Copyright 1986 Michael L. Connolly),
as distributed in the CCP4 SC program (Copyright 1999-2004 Michael Lawrence).

Generates three types of surface dots for each interface atom:
    1. Convex   — contact surface on each atom's VdW shell (outward normal)
    2. Toroidal — probe rolling between two atoms (inward normal from probe centre)
    3. Concave  — re-entrant patch where probe nestles between three atoms

All normals point away from the protein bulk, matching CCP4 SC convention.

Reference:
    Connolly, M.L. (1983) Science 221:709-713
    Lawrence, M.C. & Colman, P.M. (1993) J. Mol. Biol. 234:946-950
"""

import math
import numpy as np

# CCP4 SC defaults (from defaults.fh)
PROBE_RADIUS    = 1.7    # Å
DOT_DENSITY     = 15.0   # dots per Å²
BURIED_FLAG     = 1
ACCESSIBLE_FLAG = 0

# ---------------------------------------------------------------------------
# Radii table — exact values from CCP4 sc_radii.lib
#
# Format mirrors the Fortran match() logic:
#   - residue "***" = wildcard, matches any residue
#   - atom name ending in "*" = prefix match (e.g. "H*" matches "HA", "HB", etc.)
#   - later entries overwrite earlier wildcard matches (priority scheme)
#
# Stored as list of (residue, atom_pattern, radius) in file order so that
# specific entries correctly override wildcards.
# ---------------------------------------------------------------------------

_SC_RADII_TABLE = [
    # Wildcards first (lowest priority)
    ("***", "H",    0.50),
    ("***", "H*",   0.50),
    ("***", "H**",  0.50),
    ("***", "H***", 0.50),
    ("***", "CA",   1.85),
    ("***", "C",    1.80),
    ("***", "O",    1.60),
    ("***", "N",    1.65),
    ("***", "CB",   1.90),
    ("***", "OT*",  1.60),
    ("***", "OXT",  1.60),
    ("***", "S*",   1.90),
    ("***", "P",    1.80),
    # Residue-specific overrides (higher priority)
    ("ALA", "CB",   1.95),
    ("ARG", "NH*",  1.70),
    ("ARG", "CZ",   1.80),
    ("ARG", "NE",   1.65),
    ("ARG", "CD",   1.90),
    ("ARG", "CG",   1.90),
    ("ASN", "ND2",  1.70),
    ("ASN", "OD1",  1.60),
    ("ASN", "CG",   1.80),
    ("ASP", "OD*",  1.60),
    ("ASP", "CG",   1.80),
    ("GLN", "NE2",  1.70),
    ("GLN", "OE1",  1.60),
    ("GLN", "CD",   1.80),
    ("GLN", "CG",   1.90),
    ("GLU", "OE*",  1.60),
    ("GLU", "CD",   1.80),
    ("GLU", "CG",   1.90),
    ("GLY", "CA",   1.90),
    ("HIS", "CD2",  1.90),
    ("HIS", "NE2",  1.65),
    ("HIS", "CE1",  1.90),
    ("HIS", "ND1",  1.65),
    ("HIS", "CG",   1.80),
    ("HOH", "O**",  1.70),
    ("ILE", "CD1",  1.95),
    ("ILE", "CG1",  1.90),
    ("ILE", "CB",   1.85),
    ("ILE", "CG2",  1.95),
    ("LEU", "CD*",  1.95),
    ("LEU", "CG",   1.85),
    ("LYS", "NZ",   1.75),
    ("LYS", "CE",   1.90),
    ("LYS", "CD",   1.90),
    ("LYS", "CG",   1.90),
    ("MET", "CE",   1.95),
    ("MET", "CG",   1.90),
    ("PHE", "CD*",  1.90),
    ("PHE", "CE*",  1.90),
    ("PHE", "CZ",   1.90),
    ("PHE", "CG",   1.80),
    ("PRO", "CD",   1.90),
    ("PRO", "CG",   1.90),
    ("SER", "OG",   1.70),
    ("SUL", "S",    1.90),
    ("SUL", "O***", 1.65),
    ("THR", "CG2",  1.95),
    ("THR", "OG1",  1.70),
    ("THR", "CB",   1.85),
    ("TRP", "CE2",  1.80),
    ("TRP", "CE3",  1.90),
    ("TRP", "CD1",  1.90),
    ("TRP", "CD2",  1.80),
    ("TRP", "CZ*",  1.90),
    ("TRP", "CH2",  1.90),
    ("TRP", "NE1",  1.65),
    ("TRP", "CG",   1.80),
    ("TYR", "OH",   1.70),
    ("TYR", "CD*",  1.90),
    ("TYR", "CE*",  1.90),
    ("TYR", "CZ",   1.80),
    ("TYR", "CG",   1.80),
    ("VAL", "CG*",  1.95),
    ("VAL", "CB",   1.85),
    ("WAT", "O",    1.70),
    ("WAT", "O*",   1.70),
]

# Pre-build a lookup cache: (residue, atom) -> radius
# Uses the same priority logic as the Fortran: later matches overwrite earlier ones.
# Wildcard residue "***" matches anything; atom pattern ending in "*" is a prefix match.

def _build_radii_cache():
    """Build a (residue, atom) -> radius lookup applying Fortran priority rules."""
    cache = {}
    return cache  # populated lazily in get_radius()

_RADII_CACHE = {}

def _sc_match_atom(atom_name, pattern):
    """
    Match atom_name against pattern, where trailing '*' is a prefix wildcard.
    Mirrors the Fortran match() function.
    """
    atom = atom_name.strip().upper()
    pat  = pattern.strip().upper()
    if pat == "*":
        return True
    star = pat.find("*")
    if star == -1:
        return atom == pat
    # prefix match up to the first '*'
    return atom.startswith(pat[:star])


def get_radius(residue_name, atom_name):
    """
    Look up the CCP4 sc_radii radius for a given residue/atom pair.
    Applies the same priority scheme as the Fortran: iterate through the table
    in order; later matches overwrite earlier ones (so specific residue entries
    correctly override wildcards).

    :param residue_name: str, PDB residue name (3-letter code, e.g. 'ALA')
    :param atom_name:    str, PDB atom name (e.g. 'CA', 'OD1')
    :return: float, radius in Å
    """
    res  = residue_name.strip().upper()
    atom = atom_name.strip().upper()
    key  = (res, atom)
    if key in _RADII_CACHE:
        return _RADII_CACHE[key]

    radius = 1.80  # fallback default
    for (res_pat, atom_pat, r) in _SC_RADII_TABLE:
        res_match  = (res_pat == "***") or (res_pat.upper() == res)
        atom_match = _sc_match_atom(atom, atom_pat)
        if res_match and atom_match:
            radius = r  # later entries overwrite — same as Fortran

    _RADII_CACHE[key] = radius
    return radius


# ---------------------------------------------------------------------------
# Vector helpers (mirrors the Fortran elementary subroutines)
# ---------------------------------------------------------------------------

def _dot(a, b):      return float(np.dot(a, b))
def _cross(a, b):    return np.cross(a, b)
def _dis2(a, b):
    d = a - b; return float(np.dot(d, d))
def _dis(a, b):      return math.sqrt(max(0.0, _dis2(a, b)))
def _norm(a):
    n = math.sqrt(float(np.dot(a, a)))
    if n <= 0.0: raise ValueError("zero vector in _norm")
    return a / n
def _disptl(cen, axis, pnt):
    vec = pnt - cen; dt = _dot(vec, axis)
    return math.sqrt(max(0.0, float(np.dot(vec, vec)) - dt*dt))


# ---------------------------------------------------------------------------
# Subdivision helpers (mirrors subdiv / subarc / subcir)
# ---------------------------------------------------------------------------

def _subdiv(cen, rad, x, y, angle, density):
    if density <= 0 or rad <= 0 or angle <= 0: return [], 0.0
    delta = 1.0 / (math.sqrt(density) * rad)
    a = -delta / 2.0; pnts = []
    while True:
        a += delta
        if a > angle: break
        pnts.append(cen + rad*math.cos(a)*x + rad*math.sin(a)*y)
    if not pnts: return [], 0.0
    return pnts, rad * angle / len(pnts)

def _subarc(cen, rad, axis, density, x, v):
    y = _cross(axis, x)
    angle = math.atan2(_dot(v, y), _dot(v, x))
    if angle < 0.0: angle += 2*math.pi
    return _subdiv(cen, rad, x, y, angle, density)

def _subcir(cen, rad, axis, density):
    v1 = np.array([axis[1]**2+axis[2]**2,
                   axis[0]**2+axis[2]**2,
                   axis[0]**2+axis[1]**2], dtype=float)
    n = math.sqrt(float(np.dot(v1,v1)))
    if n > 0: v1 /= n
    if abs(_dot(v1, axis)) > 0.99: v1 = np.array([1.,0.,0.])
    v2 = _norm(_cross(axis, v1))
    x  = _norm(_cross(axis, v2))
    y  = _cross(axis, x)
    return _subdiv(cen, rad, x, y, 2.0*math.pi, density)


# ---------------------------------------------------------------------------
# Burial check
# ---------------------------------------------------------------------------

def _check_buried(pcen, atom_idx, burco, burrad_list, rp):
    """
    A dot is buried (interface-facing) if its probe centre is within
    (burrad + rp) of any atom in the opposing molecule.
    """
    for coord, r in zip(burco[atom_idx], burrad_list[atom_idx]):
        if _dis2(pcen, coord) <= (r + rp)**2:
            return True
    return False


# ---------------------------------------------------------------------------
# Main Connolly surface generator
# ---------------------------------------------------------------------------

def mds(probe_radius, atoms, radii, mol, density=DOT_DENSITY,
        residue_names=None):
    """
    Generate a Connolly molecular dot surface.

    Parameters
    ----------
    probe_radius   : float, probe radius in Å (CCP4 default: 1.7)
    atoms          : np.array (N, 3)
    radii          : np.array (N,) VdW radii
    mol            : np.array (N,) int, molecule label (1 or 2)
    density        : float, dots per Å²
    residue_names  : list of str, residue names per atom (unused, kept for
                     API compatibility — radii are passed in directly)

    Returns
    -------
    dots     : np.array (M, 3)
    normals  : np.array (M, 3)  outward unit normals
    flags    : np.array (M,) int  BURIED_FLAG or ACCESSIBLE_FLAG
    dot_mol  : np.array (M,) int  molecule label for each dot
    """
    rp     = float(probe_radius)
    natom  = len(atoms)
    atoms  = np.asarray(atoms, dtype=float)
    radii  = np.asarray(radii, dtype=float)
    mol    = np.asarray(mol, dtype=int)

    out_dots, out_normals, out_flags, out_mol = [], [], [], []

    burco       = [[] for _ in range(natom)]
    burrad_list = [[] for _ in range(natom)]
    access      = [False] * natom
    probes      = []  # (pijk, uijk, hijk, i1, i2, i3)

    # ── First loop: build neighbours, probe positions, toroidal surface ───────
    for i1 in range(natom):
        ri  = float(radii[i1])
        eri = ri + rp
        m1  = int(mol[i1])

        same_nbrs = []
        for i2 in range(natom):
            if i1 == i2: continue
            m2 = int(mol[i2])
            d2 = _dis2(atoms[i1], atoms[i2])
            bridge = ri + float(radii[i2]) + 2*rp
            if m2 == m1:
                if d2 < bridge*bridge:
                    same_nbrs.append(i2)
            else:
                if d2 < bridge*bridge:
                    burco[i1].append(atoms[i2].copy())
                    burrad_list[i1].append(float(radii[i2]))

        same_nbrs.sort(key=lambda j: _dis2(atoms[i1], atoms[j]))
        nnbr = len(same_nbrs)
        if nnbr == 0:
            access[i1] = True

        for i2 in same_nbrs:
            if i2 <= i1: continue
            rj  = float(radii[i2]); erj = rj + rp
            dij = _dis(atoms[i1], atoms[i2])
            uij = (atoms[i2] - atoms[i1]) / dij
            asymm   = (eri*eri - erj*erj) / dij
            between = abs(asymm) < dij
            tij = 0.5*(atoms[i1]+atoms[i2]) + 0.5*asymm*uij
            far = (eri+erj)**2 - dij**2
            if far <= 0: continue
            contain = dij**2 - (ri-rj)**2
            if contain <= 0: continue
            rij = 0.5*math.sqrt(far)*math.sqrt(contain)/dij

            if nnbr <= 1:
                access[i1] = access[i2] = True
            else:
                for i3 in same_nbrs:
                    if i3 <= i2: continue
                    rk = float(radii[i3]); erk = rk + rp
                    if _dis(atoms[i2], atoms[i3]) >= erj+erk: continue
                    dik = _dis(atoms[i1], atoms[i3])
                    if dik >= eri+erk: continue
                    uik = (atoms[i3] - atoms[i1]) / dik
                    dt  = max(-1., min(1., _dot(uij, uik)))
                    wijk = math.acos(dt)
                    if wijk <= 0: continue
                    swijk = math.sin(wijk)
                    if swijk <= 0: continue
                    uijk = _cross(uij, uik) / swijk
                    utb  = _cross(uijk, uij)
                    asymm_ik = (eri*eri - erk*erk) / dik
                    tik  = 0.5*(atoms[i1]+atoms[i3]) + 0.5*asymm_ik*uik
                    dt_p = _dot(uik, tik - tij)
                    bijk = tij + utb*(dt_p/swijk)
                    h2   = eri*eri - _dis2(bijk, atoms[i1])
                    if h2 <= 0:
                        dt2  = _dis2(tij, atoms[i3])
                        if dt2 < erk*erk - rij*rij: break
                        continue
                    hijk = math.sqrt(h2)
                    for isign in [+1, -1]:
                        pijk = bijk + isign*hijk*uijk
                        if any(_dis2(pijk, atoms[i4]) <= (float(radii[i4])+rp)**2
                               for i4 in same_nbrs if i4 != i2 and i4 != i3):
                            continue
                        pi1_, pi2_, pi3_ = (i1,i2,i3) if isign>0 else (i2,i1,i3)
                        probes.append((pijk.copy(), isign*uijk.copy(),
                                       hijk, pi1_, pi2_, pi3_))
                        access[i1] = access[i2] = access[i3] = True

            # Toroidal surface
            rci = rij*ri/eri; rcj = rij*rj/erj
            rb  = max(0., rij-rp)
            edens = ((rci+2*rb+rcj)/4/rij)**2 * density
            circ_pnts, ts = _subcir(tij, rij, uij, edens)
            for pij in circ_pnts:
                pij = np.array(pij)
                if any(_dis2(pij, atoms[i4]) <= (float(radii[i4])+rp)**2
                       for i4 in same_nbrs if i4 != i2):
                    continue
                access[i1] = access[i2] = True
                pi_v = (atoms[i1]-pij)/eri; pj_v = (atoms[i2]-pij)/erj
                try:  axis = _norm(_cross(pi_v, pj_v))
                except ValueError: continue
                dtq   = rp*rp - rij*rij; pcusp = dtq > 0 and between
                if pcusp:
                    dtq = math.sqrt(dtq)
                    pqi = (tij - dtq*uij - pij)/rp
                    pqj = (tij + dtq*uij - pij)/rp
                else:
                    try: pqi = _norm(pi_v + pj_v)
                    except ValueError: continue
                    pqj = pqi.copy()
                for (start_v, end_v, atom_i) in [(pi_v, pqi, i1), (pqj, pj_v, i2)]:
                    dt_ = max(-1., min(1., _dot(start_v, end_v)))
                    if abs(dt_) >= 1.: continue
                    arc, ps = _subarc(pij, rp, axis, density, start_v, end_v)
                    for pt in arc:
                        pt = np.array(pt)
                        outnml = (pij - pt) / rp
                        buried = _check_buried(pij, atom_i, burco, burrad_list, rp)
                        out_dots.append(pt); out_normals.append(outnml)
                        out_flags.append(BURIED_FLAG if buried else ACCESSIBLE_FLAG)
                        out_mol.append(m1)

    # ── Concave (re-entrant) surface ──────────────────────────────────────────
    low_idx    = [k for k, p in enumerate(probes) if p[2] < rp]
    low_coords = [probes[k][0] for k in low_idx]

    for iprb, (pijk, uijk, hijk, i1, i2, i3) in enumerate(probes):
        m1 = int(mol[i1])
        near_low = [c for k,c in zip(low_idx, low_coords)
                    if k != iprb and _dis2(pijk, c) < 4*rp*rp] if hijk < rp else []
        vp = [_norm(atoms[i] - pijk) for i in [i1, i2, i3]]
        cv = [_norm(_cross(vp[0],vp[1])),
              _norm(_cross(vp[1],vp[2])),
              _norm(_cross(vp[2],vp[0]))]
        dm, mm = -1., 0
        for m_idx in range(3):
            dt = _dot(uijk, vp[m_idx])
            if dt > dm: dm, mm = dt, m_idx
        south = -uijk
        try: axis = _norm(_cross(vp[mm], south))
        except ValueError: continue
        lats, cs = _subarc(np.zeros(3), rp, axis, density, vp[mm], south)
        for lat in lats:
            lat = np.array(lat); dt = _dot(lat, south)
            cen = dt*south; rad2 = rp*rp - dt*dt
            if rad2 <= 0: continue
            circle_pnts, ps = _subcir(cen, math.sqrt(rad2), south, density)
            for pt_raw in circle_pnts:
                pt_raw = np.array(pt_raw)
                if any(_dot(pt_raw, cv_v) >= 0 for cv_v in cv): continue
                pt = pijk + pt_raw
                if any(_dis2(pt, c) < rp*rp for c in near_low): continue
                mc, dmin_ = 0, 2*rp
                for m_idx, i in enumerate([i1, i2, i3]):
                    d = _dis(pt, atoms[i]) - float(radii[i])
                    if d < dmin_: dmin_, mc = d, i
                outnml = (pijk - pt) / rp
                buried = _check_buried(pijk, mc, burco, burrad_list, rp)
                out_dots.append(pt); out_normals.append(outnml)
                out_flags.append(BURIED_FLAG if buried else ACCESSIBLE_FLAG)
                out_mol.append(m1)

    # ── Convex (contact) surface ──────────────────────────────────────────────
    for i1 in range(natom):
        if not access[i1]: continue
        ri = float(radii[i1]); eri = ri+rp; m1 = int(mol[i1])
        same_nbrs = sorted(
            [j for j in range(natom) if j != i1 and int(mol[j]) == m1 and
             _dis2(atoms[i1], atoms[j]) < (ri+float(radii[j])+2*rp)**2],
            key=lambda j: _dis2(atoms[i1], atoms[j]))
        if not same_nbrs:
            north = np.array([0.,0.,1.]); south = np.array([0.,0.,-1.])
            eqvec = np.array([1.,0.,0.])
        else:
            i2 = same_nbrs[0]
            north = _norm(atoms[i1] - atoms[i2])
            v1 = np.array([north[1]**2+north[2]**2,
                           north[0]**2+north[2]**2,
                           north[0]**2+north[1]**2])
            n_ = math.sqrt(float(np.dot(v1,v1)))
            if n_ > 0: v1 /= n_
            if abs(_dot(v1, north)) > 0.99: v1 = np.array([1.,0.,0.])
            eqvec = _norm(_cross(north, v1))
            vql   = _cross(eqvec, north)
            rj = float(radii[i2]); erj = rj+rp
            dij = _dis(atoms[i1], atoms[i2])
            uij = (atoms[i2]-atoms[i1])/dij
            asymm = (eri*eri-erj*erj)/dij
            tij   = 0.5*(atoms[i1]+atoms[i2]) + 0.5*asymm*uij
            far   = (eri+erj)**2 - dij**2
            if far <= 0: continue
            cont  = dij**2 - (ri-rj)**2
            if cont <= 0: continue
            rij   = 0.5*math.sqrt(far)*math.sqrt(cont)/dij
            south = (tij + rij*vql - atoms[i1]) / eri
        lats, cs = _subarc(np.zeros(3), ri, eqvec, density, north, south)
        for lat in lats:
            lat = np.array(lat); dt = _dot(lat, north)
            cen = atoms[i1] + dt*north; rad2 = ri*ri - dt*dt
            if rad2 <= 0: continue
            circle_pnts, ps = _subcir(cen, math.sqrt(rad2), north, density)
            for pt_raw in circle_pnts:
                pt_raw = np.array(pt_raw)
                pcen = atoms[i1] + (eri/ri)*(pt_raw - atoms[i1])
                if any(_dis2(pcen, atoms[i4]) <= (float(radii[i4])+rp)**2
                       for i4 in same_nbrs[1:]):
                    continue
                outnml = _norm(pt_raw - atoms[i1])
                buried = _check_buried(pcen, i1, burco, burrad_list, rp)
                out_dots.append(pt_raw); out_normals.append(outnml)
                out_flags.append(BURIED_FLAG if buried else ACCESSIBLE_FLAG)
                out_mol.append(m1)

    if not out_dots:
        return (np.empty((0,3)), np.empty((0,3)),
                np.empty(0, dtype=int), np.empty(0, dtype=int))

    return (np.array(out_dots), np.array(out_normals),
            np.array(out_flags, dtype=int), np.array(out_mol, dtype=int))


# ---------------------------------------------------------------------------
# Trim band (mirrors Fortran trim subroutine exactly)
# ---------------------------------------------------------------------------

def trim(dots, flags, band):
    """
    Remove buried dots within `band` Å of any accessible dot on the same surface.
    Returns boolean mask — True = dot survives.
    """
    from scipy.spatial import cKDTree
    buried_mask = flags == BURIED_FLAG
    acc_mask    = flags == ACCESSIBLE_FLAG
    if not acc_mask.any():
        return buried_mask.copy()
    acc_tree = cKDTree(dots[acc_mask])
    keep = np.zeros(len(dots), dtype=bool)
    for i in np.where(buried_mask)[0]:
        if len(acc_tree.query_ball_point(dots[i], band)) == 0:
            keep[i] = True
    return keep
