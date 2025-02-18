"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid typically contains a corrin macrocycle with cobalt.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cobalt atom
    cobalt_atom = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'Co':
            cobalt_atom = atom
            break
    if not cobalt_atom:
        return False, "Cobalt atom not found"
    
    # Define a SMARTS pattern for a corrin-like macrocyclic ring, simplified as complex structure
    # Note: This is a rough representation and might need further refinement
    corrin_pattern = Chem.MolFromSmarts("[#6;R][#7;R][$([#7]=[#6]-[#6])]-[$([#6]=[#6]-[#7])][#6;R]")
    
    # Check if the molecule contains this pattern
    if not mol.HasSubstructMatch(corrin_pattern):
        return False, "Corrin-like macrocycle pattern not found"

    return True, "Contains cobalt and corrin-like macrocycle consistent with corrinoid structure"