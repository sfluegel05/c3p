"""
Classifies: CHEBI:33521 metal atom
"""
from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a metal atom, False otherwise
        str: Reason for classification
    """
    # List of metal elements by symbol
    metal_symbols = [
        "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti",
        "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb",
        "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
        "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
        "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
        "Pb", "Bi", "Po", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np",
        "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Fl",
        "Lv"
    ]

    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure it represents a single atom
    if mol.GetNumAtoms() != 1:
        return False, "SMILES string does not represent a single atom"

    # Extract the element symbol
    atom = mol.GetAtomWithIdx(0)
    element_symbol = atom.GetSymbol()

    # Check if the element is in the list of metal symbols
    if element_symbol in metal_symbols:
        return True, f"{element_symbol} is a metal atom"
  
    return False, f"{element_symbol} is not a metal atom"