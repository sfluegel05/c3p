"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin consists of four pyrrole rings interconnected via methine bridges,
    forming a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is identified as a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for porphyrin core structure.
    # This pattern focuses on the 18-atom macrocycle composed of four linked pyrrole units.
    porphyrin_pattern = Chem.MolFromSmarts('n1c(cc2ncc(c2)c3ncc(c3)c4ncc(c4)c1)')
    
    if porphyrin_pattern is None:
        return None, "Error in defining porphyrin pattern"

    # Check if the molecule matches the porphyrin macrocyclic structure pattern.
    if mol.HasSubstructMatch(porphyrin_pattern):
        return True, "Matches porphyrin macrocyclic structure"

    return False, "Does not match porphyrin structure"