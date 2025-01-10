"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin consists of four pyrrole rings interconnected via methine bridges.

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
    
    # Attempt to define a generalized SMARTS pattern for the core porphyrin structure:
    # The pattern attempts to match the four pyrrole rings interconnected,
    # Considering possible post coordination bonds and varying substitutions
    porphyrin_pattern = Chem.MolFromSmarts('C1=Cc2ccn([FeZnMg]*-*)c2C=c2[nH]c3ccc[nH]c3C=c1.N4')

    # Check if the molecule matches the generalized porphyrin core pattern
    if mol.HasSubstructMatch(porphyrin_pattern):
        return True, "Matches porphyrin macrocyclic structure"

    return False, "Does not match porphyrin structure"