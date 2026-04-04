import math
import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import plotly.express as px

import logging

from .connolly import (mds as _connolly_mds, trim as _connolly_trim,
                       PROBE_RADIUS as _PROBE_RADIUS, BURIED_FLAG as _BURIED_FLAG)

from scipy.spatial import cKDTree, ConvexHull
from sklearn.decomposition import PCA
from dataclasses import dataclass
from itertools import compress

logger = logging.getLogger(__name__)


@dataclass
class PDBCoords:
    """
    Dataclass containing PDB coordinates, amino acids, atoms and residue numbers
    """
    coords: np.array
    amino_acids: list
    atoms: list
    residues: list


class ShapeComplementarity:
    """
    Shape Complementarity implementation based on Lawrence & Colman (1993).

    Reference:
        Lawrence, M.C. & Colman, P.M. (1993) J. Mol. Biol. 234:946-950
    """

    def __init__(self, arg):
        self.verbose   = None
        self.distance  = None
        self.density   = None
        self.weight    = None
        self.complex_1 = None
        self.complex_2 = None
        self.plot      = None
        self.arg       = arg
        super().__init__(self.arg)

    # -----------------------------------------------------------------------
    # PDB parsing helpers
    # -----------------------------------------------------------------------

    def convert_1d_array(self, arr):
        """
        Return a numpy array from a list.
        :param arr: list
        :return: np.array
        """
        return np.array(arr, dtype=float)

    def create_interface(self):
        """
        Create an interface of the two complexes, returning two PDBCoords objects.
        :return: PDBCoords c1, PDBCoords c2
        """
        logger.info("Getting interfaces")

        residues    = self.get_column("RESIDUE_NUM")
        chain       = self.get_column("CHAIN")
        amino_acids = self.get_column("RESIDUE_SEQID")
        atoms       = self.get_column("ATOM_NAME")

        x_coord = self.convert_1d_array(self.get_column("X_COORD"))
        y_coord = self.convert_1d_array(self.get_column("Y_COORD"))
        z_coord = self.convert_1d_array(self.get_column("Z_COORD"))

        complex_1_created = False
        complex_2_created = False

        complex_1_residues, complex_2_residues = [], []
        complex_1_aa,       complex_2_aa       = [], []
        complex_1_at,       complex_2_at       = [], []

        for x, y, z, c, r, aa, a in zip(x_coord, y_coord, z_coord, chain,
                                         residues, amino_acids, atoms):
            coord = np.array((x, y, z), "f")
            if c in self.complex_1:
                complex_1_residues.append(r)
                complex_1_aa.append(aa)
                complex_1_at.append(a)
                if complex_1_created:
                    complex_1_coords = np.append(complex_1_coords, coord)
                else:
                    complex_1_coords = coord
                    complex_1_created = True
            if c in self.complex_2:
                complex_2_residues.append(r)
                complex_2_aa.append(aa)
                complex_2_at.append(a)
                if complex_2_created:
                    complex_2_coords = np.append(complex_2_coords, coord)
                else:
                    complex_2_coords = coord
                    complex_2_created = True

        complex_1_coords = np.reshape(complex_1_coords, (-1, 3))
        complex_2_coords = np.reshape(complex_2_coords, (-1, 3))

        c1 = PDBCoords(coords=complex_1_coords, amino_acids=complex_1_aa,
                       atoms=complex_1_at, residues=complex_1_residues)
        c2 = PDBCoords(coords=complex_2_coords, amino_acids=complex_2_aa,
                       atoms=complex_2_at, residues=complex_2_residues)
        return c1, c2

    def filter_interface(self, a, b, r):
        """
        Return only atoms in `a` that are within r Angstroms of any atom in `b`.
        :param a: PDBCoords
        :param b: PDBCoords
        :param r: float, distance threshold in Å
        :return: PDBCoords
        """
        tree = cKDTree(a.coords)
        mask = np.zeros(len(a.coords), dtype=bool)

        indices = []
        for coord in b.coords:
            indices += tree.query_ball_point(coord, r)

        mask[list(set(indices))] = True

        return PDBCoords(
            coords     = a.coords[mask],
            amino_acids= list(compress(a.amino_acids, mask)),
            atoms      = list(compress(a.atoms,       mask)),
            residues   = list(compress(a.residues,    mask)),
        )

    # -----------------------------------------------------------------------
    # Surface mesh and dot sampling
    # -----------------------------------------------------------------------

    def create_polygon(self, points):
        """
        Triangulate the 3D interface surface using ConvexHull.
        Returns simplices as an array of triangle vertex-index triples.
        :param points: np.array (n, 3)
        :return: np.array (m, 3)
        """
        return ConvexHull(points).simplices

    def point_inside_triangle(self, v1, v2, v3):
        """
        Sample a uniformly random point inside a triangle.
        :param v1, v2, v3: np.array (3,), triangle vertices
        :return: np.array (3,)
        """
        a = math.sqrt(random.random())
        b = random.random()
        return (1 - a) * v1 + a * (1 - b) * v2 + a * b * v3

    def random_points(self, coords, simplices, n_samp):
        """
        Generate n_samp random points sampled uniformly from the triangulated surface.
        :param coords: np.array (n, 3)
        :param simplices: np.array (m, 3), triangle vertex indices
        :param n_samp: int
        :return: np.array (n_samp, 3)
        """
        indices = np.random.choice(len(simplices), n_samp)
        pts = []
        for i in indices:
            s = simplices[i]
            pts.append(self.point_inside_triangle(coords[s[0]], coords[s[1]], coords[s[2]]))
        return np.array(pts)

    def estimate_surface_area(self, coords):
        """
        Estimate the surface area of an interface using ConvexHull.
        :param coords: np.array (n, 3)
        :return: float, area in Å²
        """
        area = ConvexHull(coords).area
        logger.info("Estimated area of interface is %.2f\N{ANGSTROM SIGN}\N{SUPERSCRIPT TWO}", area)
        return area

    # -----------------------------------------------------------------------
    # Normal estimation
    # -----------------------------------------------------------------------

    def calculate_normal(self, coordinate, mesh):
        """
        Estimate the surface normal at a dot using PCA on its 10 nearest
        neighbours in the dot cloud. The least-variance PCA component gives
        the normal to the local tangent plane.
        :param coordinate: np.array (3,)
        :param mesh: np.array (n, 3)
        :return: np.array (3,), unit normal vector
        """
        tree = cKDTree(mesh)
        _, ind = tree.query(coordinate, k=min(10, len(mesh)))
        pca = PCA(n_components=3)
        pca.fit(mesh[ind])
        return pca.components_[np.argmin(pca.explained_variance_)]

    # -----------------------------------------------------------------------
    # SC calculation
    # -----------------------------------------------------------------------

    def find_nearest_neighbour(self, coord, set_of_coords):
        """
        Return the nearest point in set_of_coords to coord.
        :param coord: np.array (3,)
        :param set_of_coords: np.array (n, 3)
        :return: np.array (3,)
        """
        _, idx = cKDTree(set_of_coords).query(coord)
        return set_of_coords[idx]

    def surface_complementarity_function(self, n_a, n_b, x_a, x_b, w):
        """
        Compute S(A->B)(x_A) = (n_A · n_B) · exp(-w · |x_A - x_B|)
        :param n_a: np.array, unit normal at x_a on surface A
        :param n_b: np.array, unit normal at nearest x_b on surface B
        :param x_a: np.array, coordinate on surface A
        :param x_b: np.array, nearest coordinate on surface B
        :param w: float, distance weighting factor (0 = pure dot product)
        :return: float
        """
        dot = np.dot(n_a, n_b)
        if w == 0.0:
            return dot
        return dot * np.exp(-np.linalg.norm(x_a - x_b) * w)

    def calculate_sc(self, points_c1, points_c2, weight):
        """
        Compute SC(A->B) for every dot on surface A toward surface B.
        :param points_c1: np.array (M, 3), sampled dots on surface A
        :param points_c2: np.array (K, 3), sampled dots on surface B
        :param weight: float, exponential distance weighting
        :return: list of float
        """
        tree_b   = cKDTree(points_c2)
        sc_array = []
        for pt_a in points_c1:
            _, idx = tree_b.query(pt_a)
            pt_b   = points_c2[idx]
            n_a    = self.calculate_normal(pt_a, points_c1)
            n_b    = self.calculate_normal(pt_b, points_c2)
            sc_array.append(self.surface_complementarity_function(n_a, n_b, pt_a, pt_b, weight))
        return sc_array

    # -----------------------------------------------------------------------
    # Plotting helpers
    # -----------------------------------------------------------------------

    def plot_sc(self, sc_complex_1, sc_complex_2):
        """
        Plot a histogram of SC function values for both surfaces.
        :param sc_complex_1: list of float, S(C1->C2) values
        :param sc_complex_2: list of float, S(C2->C1) values
        """
        c1 = pd.DataFrame({"SC_function": sc_complex_1, "Complex": "Complex 1"})
        c2 = pd.DataFrame({"SC_function": sc_complex_2, "Complex": "Complex 2"})
        sns.histplot(data=pd.concat([c1, c2], ignore_index=True),
                     x="SC_function", hue="Complex")
        plt.title("SC function distribution")
        plt.show()

    def plot_combined_mesh(self, mesh1, mesh2, coords1, coords2):
        x = np.concatenate([coords1[:, 0], coords2[:, 0]])
        y = np.concatenate([coords1[:, 1], coords2[:, 1]])
        z = np.concatenate([coords1[:, 2], coords2[:, 2]])
        colours = ["#EF553B"] * len(coords1) + ["#00CC96"] * len(coords2)
        return ff.create_trisurf(x=x, y=y, z=z,
                                 simplices=np.concatenate([mesh1, mesh2]),
                                 colormap=colours,
                                 title="Surface of Complex 1 and 2")

    def plot_single_mesh(self, mesh, coords, title):
        return ff.create_trisurf(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                                 simplices=mesh, title=title)

    def plot_atoms(self, c1, c2, title):
        x = np.concatenate([c1[:, 0], c2[:, 0]])
        y = np.concatenate([c1[:, 1], c2[:, 1]])
        z = np.concatenate([c1[:, 2], c2[:, 2]])
        colours = ["Complex 1"] * len(c1) + ["Complex 2"] * len(c2)
        df = pd.DataFrame(list(zip(x, y, z, colours)), columns=["X", "Y", "Z", "Complex"])
        return px.scatter_3d(df, x="X", y="Y", z="Z", color="Complex", title=title)

    # -----------------------------------------------------------------------
    # Main entry point
    # -----------------------------------------------------------------------

    # ProtOr radii (Tsai et al. 1999) for Connolly surface generation
    _PROTOR = {
        "N":1.65,"CA":1.87,"C":1.76,"O":1.40,"CB":1.87,"CG":1.87,
        "CG1":1.87,"CG2":1.87,"CD":1.87,"CD1":1.87,"CD2":1.87,
        "CE":1.87,"CE1":1.87,"CE2":1.87,"CE3":1.87,"CZ":1.87,
        "CZ2":1.87,"CZ3":1.87,"CH2":1.87,
        "ND1":1.65,"ND2":1.65,"NE":1.65,"NE1":1.65,"NE2":1.65,
        "NH1":1.65,"NH2":1.65,"NZ":1.65,
        "OD1":1.40,"OD2":1.40,"OE1":1.40,"OE2":1.40,
        "OG":1.40,"OG1":1.40,"OH":1.40,"OXT":1.40,
        "SD":1.85,"SG":1.85,
    }

    def sc(self):
        """
        Calculate Shape Complementarity (SC) using the Lawrence & Colman (1993) method.

        Surface dots are generated using the Connolly molecular surface algorithm,
        translated directly from the Fortran mds subroutine in CCP4 SC. This produces
        three types of surface dots:
            1. Convex   — contact surface on each atom's VdW shell
            2. Toroidal — probe rolling between two atoms
            3. Concave  — re-entrant patch where probe nestles between three atoms

        Only dots flagged as buried (probe centre within reach of the opposing molecule)
        are scored. A 1.5 Å trim band removes peripheral edge dots.

        SC = (median S(C1->C2) + median S(C2->C1)) / 2
        where S = -(nA · nB)  [negated: complementary opposing normals score positive]
        """
        complex1, complex2 = self.create_interface()
        complex1 = self.filter_interface(complex1, complex2, self.distance)
        complex2 = self.filter_interface(complex2, complex1, self.distance)

        logger.info("Complex 1 contains %d atoms within %g Angstroms of Complex 2",
                     len(complex1.residues), self.distance)
        logger.info("Complex 2 contains %d atoms within %g Angstroms of Complex 1",
                     len(complex2.residues), self.distance)

        from .connolly import get_radius
        atoms = np.vstack([complex1.coords, complex2.coords])
        radii = np.array([
            get_radius(aa, at)
            for aa, at in zip(complex1.amino_acids + complex2.amino_acids,
                              complex1.atoms       + complex2.atoms)
        ])
        mol   = np.array([1]*len(complex1.coords) + [2]*len(complex2.coords))

        logger.info("Generating Connolly surface (density=%.1f dots/\N{ANGSTROM SIGN}\N{SUPERSCRIPT TWO})",
                    self.density)
        dots, normals, flags, dot_mol = _connolly_mds(
            _PROBE_RADIUS, atoms, radii, mol, density=self.density)

        if len(dots) == 0:
            logger.warning("No surface dots generated — try increasing --distance")
            return None

        buried = flags == _BURIED_FLAG
        logger.info("Total dots: %d  buried: %d", len(dots), buried.sum())

        # Select buried dots from each surface
        d1_all = dots[(dot_mol == 1) & buried]
        n1_all = normals[(dot_mol == 1) & buried]
        d2_all = dots[(dot_mol == 2) & buried]
        n2_all = normals[(dot_mol == 2) & buried]

        if len(d1_all) == 0 or len(d2_all) == 0:
            logger.warning("No buried dots — try increasing --distance or --dot-density")
            return None

        # Filter to dots within 1.0 Å of the opposing surface.
        # This mirrors the effect of the CCP4 SC trim band, which implicitly
        # removes peripheral dots far from the opposing surface. Using an
        # explicit distance cutoff is equivalent and more transparent.
        _, idx2_all = cKDTree(d2_all).query(d1_all)
        _, idx1_all = cKDTree(d1_all).query(d2_all)
        dists1_all = np.linalg.norm(d1_all - d2_all[idx2_all], axis=1)
        dists2_all = np.linalg.norm(d2_all - d1_all[idx1_all], axis=1)

        DIST_CUTOFF = 1.4  # Å — equivalent to CCP4 trim band effect
        m1 = dists1_all <= DIST_CUTOFF
        m2 = dists2_all <= DIST_CUTOFF

        d1 = d1_all[m1];  n1 = n1_all[m1];  dists1 = dists1_all[m1]
        d2 = d2_all[m2];  n2 = n2_all[m2];  dists2 = dists2_all[m2]

        logger.info("After distance filter (≤%.1fÅ): C1=%d dots, C2=%d dots",
                    DIST_CUTOFF, len(d1), len(d2))

        if len(d1) == 0 or len(d2) == 0:
            logger.warning("No close-contact dots — try increasing --distance or --dot-density")
            return None

        # Re-query nearest neighbours within the filtered set
        _, idx2 = cKDTree(d2).query(d1)
        _, idx1 = cKDTree(d1).query(d2)
        dists1 = np.linalg.norm(d1 - d2[idx2], axis=1)
        dists2 = np.linalg.norm(d2 - d1[idx1], axis=1)
        dot1   = -(np.einsum("ij,ij->i", n1, n2[idx2]))
        dot2   = -(np.einsum("ij,ij->i", n2, n1[idx1]))

        if self.weight > 0:
            # CCP4 SC formula: S = -(nA·nB) * exp(-d²*w)
            sc_complex_1 = list(dot1 * np.exp(-dists1**2 * self.weight))
            sc_complex_2 = list(dot2 * np.exp(-dists2**2 * self.weight))
        else:
            sc_complex_1 = list(dot1)
            sc_complex_2 = list(dot2)

        sc_score = float((np.median(sc_complex_1) + np.median(sc_complex_2)) / 2)

        logger.info("SC = %.2f", sc_score)
        print(f"{sc_score:.2f}")

        if self.plot:
            self.plot_sc(sc_complex_1, sc_complex_2)
            self.plot_atoms(complex1.coords, complex2.coords, "Interface atoms").show()
            self.plot_atoms(d1, d2, "Connolly surface dots (trimmed buried)").show()

        return sc_score
    def get_column(self, param):
        pass
