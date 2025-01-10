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
    
    # Comprehensive list of metallic elements by their symbols
    metal_symbols = set([
        'Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V',  
        'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 
        'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 
        'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 
        'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 
        'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 
        'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Og', 'Fr', 
        'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 
        'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
    ])

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the SMILES represents a single atom irrespective of charge or isotopes
    if mol.GetNumAtoms() != 1:
        return False, "SMILES does not represent a single atom"
    
    # Get the atom, verify its elemental symbol
    atom = mol.GetAtomWithIdx(0)
    symbol = atom.GetSymbol()

    # Verify if the atom symbol is in the list of metallic elements
    if symbol in metal_symbols:
        # Clarify using isotope information where applicable
        isotope = atom.GetIsotope()
        isotope_info = f" (isotope: {isotope})" if isotope else ""
        return True, f"{symbol}{isotope_info} is a metal atom"
    else:
        return False, f"{symbol} is not a metal atom"