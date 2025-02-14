"""
Classifies: CHEBI:33521 metal atom
"""
from rdkit import Chem
from rdkit.Chem import PeriodicTable

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    A metal atom is defined as a single atom of a metallic element.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): (True, reason) if molecule is a metal atom, (False, reason) otherwise.
                          If parsing fails return (None, None)
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None

        # Check if the molecule is a single atom
        if mol.GetNumAtoms() != 1:
            return False, "Not a single atom."

        atom = mol.GetAtomWithIdx(0)
        symbol = atom.GetSymbol()

        # Get atomic number from the symbol using PeriodicTable
        pt = PeriodicTable.GetPeriodicTable()
        atomic_num = pt.GetAtomicNumber(symbol)


        # Check if the atom is a metal
        if symbol in ["Li","Na","K","Rb","Cs","Fr", "Be","Mg","Ca","Sr","Ba","Ra",
                    "Sc","Y","La","Ac",
                    "Ti","Zr","Hf","Rf",
                    "V","Nb","Ta","Db",
                    "Cr","Mo","W","Sg",
                    "Mn","Tc","Re","Bh",
                    "Fe","Ru","Os","Hs",
                    "Co","Rh","Ir","Mt",
                    "Ni","Pd","Pt","Ds",
                    "Cu","Ag","Au","Rg",
                    "Zn","Cd","Hg","Cn",
                    "Al","Ga","In","Tl",
                     "Ge","Sn","Pb",
                    "Sb","Bi",
                    "Po"]:
            return True, "Single metal atom"
        elif atomic_num > 0 :
            return False, f"Single atom, but not a metal: {symbol}"
        else:
             return False, "Invalid atom"
    except Exception as e:
        return None, None