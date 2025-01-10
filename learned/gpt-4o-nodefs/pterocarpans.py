"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans have a polycyclic structure often described
    as a dihydrobenzofuro[3,2-c]chromene core.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True and reason if the molecule is a pterocarpan, False and reason otherwise
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern that captures the core dihydrobenzofurochromene structure
    # The pattern focuses on the typical 6a,11a-dihydrobenzofuro[3,2-c]chromene core
    # Adjust the pattern based on further studies or example compounds
    pterocarpan_pattern = Chem.MolFromSmarts("C1OC2C(C3=CC=CC=C3C(O4)C21)C4")
    
    if pterocarpan_pattern is None:
        return (None, "Error in constructing SMARTS pattern")

    # Check for the presence of the pterocarpan core structure within the molecule
    if mol.HasSubstructMatch(pterocarpan_pattern):
        return True, "Contains the pterocarpan core structure"
    
    return False, "Does not contain the pterocarpan core structure"