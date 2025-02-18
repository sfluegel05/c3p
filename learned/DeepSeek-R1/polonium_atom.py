"""
Classifies: CHEBI:33313 polonium atom
"""
"""
Classifies: CHEBI:33317 polonium atom
"""
from rdkit import Chem

def is_polonium_atom(smiles: str):
    """
    Determines if a SMILES string represents a single polonium atom (any isotope).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if it's a single polonium atom, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    
    # Check if SMILES is valid
    if not mol:
        return False, "Invalid SMILES"
    
    # Check for single atom
    if mol.GetNumAtoms() != 1:
        return False, "Not a single atom"
    
    # Get the only atom in the molecule
    atom = mol.GetAtomWithIdx(0)
    
    # Check if it's polonium (atomic number 84)
    if atom.GetAtomicNum() != 84:
        return False, "Atom is not polonium"
    
    return True, "Single polonium atom"