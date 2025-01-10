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
        bool: True if molecule is a metal atom, False otherwise
        str: Reason for classification
    """
    
    # List of metallic elements by their symbols
    metal_symbols = set([
        'Li', 'Be', 'Na', 'Mg', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 
        'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 
        'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 
        'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
        'Po', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 
        'Fl', 'Mc', 'Lv', 'Og', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np',
        'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
    ])

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the SMILES represents a single atom despite charge or isotopes
    if mol.GetNumAtoms() != 1:
        return False, "SMILES does not represent a single atom"
    
    # Get the atom
    atom = mol.GetAtomWithIdx(0)
    symbol = atom.GetSymbol()

    # Check isotopes and charge
    isotope = atom.GetIsotope()
    charge = atom.GetFormalCharge()

    # Verify if the atom symbol is in the list of metallic elements
    if symbol in metal_symbols:
        # Report considering isotopes and charge
        isotope_info = f" (isotope: {isotope})" if isotope else ""
        charge_info = f" (charge: {charge:+})" if charge != 0 else ""
        return True, f"{symbol}{isotope_info}{charge_info} is a metal atom"
    else:
        return False, f"{symbol} is not a metal atom"