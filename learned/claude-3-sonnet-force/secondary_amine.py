"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:35550 secondary amine
A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for secondary amine
    pattern = Chem.MolFromSmarts("[NX3H2]([CX4,cX3])[CX4,cX3]")
    
    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains a secondary aliphatic amine group"
    else:
        return False, "Does not contain a secondary aliphatic amine group"