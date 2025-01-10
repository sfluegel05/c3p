"""
Classifies: CHEBI:33313 polonium atom
"""
from rdkit import Chem

def is_polonium_atom(smiles: str):
    """
    Determines if a molecule is a polonium atom based on its SMILES string.
    A polonium atom is denoted by a SMILES string of the form [massPo], where
    "mass" is the atomic mass of the polonium isotope.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polonium atom, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule contains exactly one heavy atom
    num_atoms = mol.GetNumAtoms()
    if num_atoms != 1:
        return False, "SMILES does not represent a single atom"
    
    # Retrieve the atom and check its element and isotope
    atom = mol.GetAtomWithIdx(0)
    if atom.GetSymbol() != "Po":
        return False, "The atom is not polonium"

    # Check if isotope information is present
    isotope = atom.GetIsotope()
    if isotope == 0:
        return False, "No isotope information for polonium"

    return True, f"Polonium-atom with isotope {isotope}"