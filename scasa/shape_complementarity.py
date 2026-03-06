import math
import random

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import plotly.express as px

from scipy.spatial import cKDTree, ConvexHull
from sklearn.decomposition import PCA
from dataclasses import dataclass
from itertools import compress


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
    Shape Complementarity implementation based on Lawrence & Colman (1993)
    but with some python library shortcuts

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
        if self.verbose:
            print("Getting interfaces")

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

        return PDBCoords(coords = a.coords[mask], amino_acids = list(compress(a.amino_acids, mask)),
                         atoms = list(compress(a.atoms, mask)), residues = list(compress(a.residues, mask)))


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
        if self.verbose:
            print(f"Estimated area of complex 1's face is {area:.2f}\N{ANGSTROM SIGN}\N{SUPERSCRIPT TWO}")
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
        c1 = pd.DataFrame(sc_complex_1, columns=["SC_function"])
        c1["Complex"] = "Complex 1"
        c2 = pd.DataFrame(sc_complex_2, columns=["SC_function"])
        c2["Complex"] = "Complex 2"
        sns.histplot(data=pd.concat([c1, c2]), x="SC_function", hue="Complex")
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

    def sc(self):
        """
        Calculate Shape Complementarity (SC) using the Lawrence & Colman (1993) method.

        Interface atoms are filtered by --distance. A ConvexHull triangulation of
        each interface is randomly sampled to produce surface dots. Per-dot normals
        are estimated via PCA on the 10 nearest dot neighbours.

        SC = (median S(C1->C2) + median S(C2->C1)) / 2
        """
        complex1, complex2 = self.create_interface()
        complex1 = self.filter_interface(complex1, complex2, self.distance)
        complex2 = self.filter_interface(complex2, complex1, self.distance)

        if self.verbose:
            print(f"Complex 1 contains {len(complex1.residues)} atoms "
                  f"within {self.distance} Angstroms of Complex 2")
            print(f"Complex 2 contains {len(complex2.residues)} atoms "
                  f"within {self.distance} Angstroms of Complex 1")

        area_1 = self.estimate_surface_area(complex1.coords)
        area_2 = self.estimate_surface_area(complex2.coords)

        simplices_c1 = self.create_polygon(complex1.coords)
        simplices_c2 = self.create_polygon(complex2.coords)

        n_dots_1 = round(self.density * area_1)
        n_dots_2 = round(self.density * area_2)

        points_c1 = self.random_points(complex1.coords, simplices_c1, n_dots_1)
        points_c2 = self.random_points(complex2.coords, simplices_c2, n_dots_2)

        if self.verbose:
            print("Calculating SC for both complexes")

        sc_complex_1 = self.calculate_sc(points_c1, points_c2, self.weight)
        sc_complex_2 = self.calculate_sc(points_c2, points_c1, self.weight)

        sc_score = (np.median(sc_complex_1) + np.median(sc_complex_2)) / 2

        if self.verbose:
            print(f"SC = {sc_score:.2f}")
        else:
            print(f"{sc_score:.2f}")

        return sc_score

    def get_column(self, param):
        pass
