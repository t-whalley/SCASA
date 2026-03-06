from .available_surface_area import SurfaceArea
from .shape_complementarity import ShapeComplementarity

from pathlib import Path


class ChainNotFoundException(Exception):
    """
    Custom exception raised when there are no variants detected
    """

    def __init__(self, arg):
        self.arg = arg
        super().__init__(self.arg)


class Complex(SurfaceArea, ShapeComplementarity):
    """
    Object for defining and holding a PDB complex
    PDBFile: path to PDB file
    Complex1: Complex of chains, supplied as a string
    Complex2: As Complex1, but if left as None, it will be assumed to be all remaining chains in PDBFile - Complex1
    Verbose: Extra logging
    """

    def __init__(self, pdb_file, complex_1, complex_2=None, verbose=False, tmp_directory="/tmp", distance=8,
                 density=1.5, weight=0.0, plot=False):
        super().__init__()
        self.pdb_ranges = {"ATOM": range(0, 4),
                           "SERIAL": range(6, 11),
                           "ATOM_NAME": range(12, 16),
                           "ALT_LOC": range(16, 17),
                           "RESIDUE_SEQID": range(17, 20),
                           "CHAIN": range(21, 22),
                           "RESIDUE_NUM": range(22, 26),
                           "INSERTION": range(26, 27),
                           "X_COORD": range(30, 38),
                           "Y_COORD": range(38, 46),
                           "Z_COORD": range(46, 54),
                           "OCCUPANCY": range(54, 60),
                           "TEMPERATURE": range(60, 66),
                           "SEGMENT": range(72, 76),
                           "ELEMENT": range(76, 78)}

        self.pdb_file = pdb_file
        self.complex_1 = complex_1
        self.complex_2 = complex_2
        self.chains = None
        self.verbose = verbose
        self.tmp_directory = tmp_directory

        # variables inherited from SC
        self.density = density
        self.distance = distance
        self.weight = weight
        self.plot = plot
        self.probe = 1.7

        # for time being assume pdb file is .pdb
        # todo change extension pattern
        self.pdb_name = pdb_file.split("/")[-1].replace(".pdb", "")

        self.verify_chains()
        self.chain_string_to_list()

        self.complex_1_asa_df = None
        self.complex_2_asa_df = None
        self.complex_1_2_asa_df = None

        if self.verbose:
            print("PDB successfully validated")
            print(f"Complex 1: {''.join(self.complex_1)}")
            print(f"Complex 2: {''.join(self.complex_2)}")

    def verify_chains(self):
        """
        Function to check PDB file exists and chains are in PDB file; then populate Complex2 if left empty
        :return: None
        """
        if not Path(self.pdb_file).is_file():
            raise FileNotFoundError(f"{self.pdb_file} does not exist.")

        self.get_all_chains()

        # Check complex 1's chains
        for chain in self.complex_1:
            if chain not in self.chains:
                raise ChainNotFoundException(f" Chain {chain} in complex_1 was not found in PDB file {self.pdb_file}")

        # and now complex 2
        if self.complex_2 is None:
            self.complex_2 = ""
            for chain in self.chains:
                if chain not in list(self.complex_1):
                    self.complex_2 += chain
            else:
                for chain in self.complex_2:
                    if chain not in self.chains:
                        raise ChainNotFoundException(
                            f" Chain {chain} in complex_2 was not found in PDB file {self.pdb_file}")

    def get_column(self, key):
        """
        Function to extract columns from PDB file based on space indexing.
        param key: key for dictionary containing deliminations of columns in POB file
        :return: list_of_col: list of rows found in that PDB column, based on indexes in pdb_ranges
        """
        list_of_col = []

        with open(self.pdb_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    indexes = list(self.pdb_ranges[key])
                    contents = "".join([line[x] for x in indexes]).strip()
                    list_of_col.append(contents)

        return list_of_col

    def get_all_chains(self):
        """
        Get all chains in a PDB file that are unique.
        :return:
        """
        chains = set(self.get_column("CHAIN"))
        self.chains = sorted(list(set(chains)))
        if self.verbose:
            print(f"PDB file contains chains: {self.chains}")

    def chain_string_to_list(self):
        """
        convert complex string to list for consistency after possibly parsing them.
        :return:
        """
        self.complex_1 = list(self.complex_1)
        self.complex_2 = list(self.complex_2)

    def subset_pdb(self, chains, pdb_out):
        outfile = open(pdb_out, "w")

        with open(self.pdb_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    indexes = list(self.pdb_ranges["CHAIN"])
                    contents = "".join([line[x] for x in indexes]).strip()
                    if contents in chains:
                        outfile.write(line)

    def create_sub_pdbs(self):
        self.subset_pdb(self.complex_1, "tmp/complex_1.pdb")
        self.subset_pdb(self.complex_2, "tmp/complex_2.pdb")
        self.subset_pdb(self.complex_1 + self.complex_2, "tmp/complex_1_2.pdb")
