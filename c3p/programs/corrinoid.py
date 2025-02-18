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
    
    # Define a more complex SMARTS pattern for detecting corrin-like macrocyclic structure
    # Focus on connected ring systems with multiple heteroatoms coordinated around a central cobalt
    corrin_pattern = Chem.MolFromSmarts("""
    [Co]1:[#7]:[#6]:[#6,=#7]-[#6](-[#7]1):[#6]=[#6]-[#6](=[#7])-[#6](-[#6](-[#7]):[Co]1:[#7]:[#6]:[#6,=#7]-[#6](-[#7]1):[#6]=[#6]-[#6]-S(=O)(=O)[O-])
    """)
    
    # Check if the molecule contains this pattern
    if not mol.HasSubstructMatch(corrin_pattern):
        return False, "Corrin-like macrocycle pattern not found"

    return True, "Contains cobalt and corrin-like macrocycle consistent with a corrinoid structure"