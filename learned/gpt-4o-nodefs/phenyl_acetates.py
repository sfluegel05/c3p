"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate derivative based on its SMILES string.
    A phenyl acetate has a phenyl group that can possibly hold an acetate group (via oxygen or similar linkage).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Updated SMARTS to identify phenyl-acetate scenario in various configurations
    # Account for substitutions on phenyl ring, and single ester linkage
    phenyl_acetate_patterns = [
        Chem.MolFromSmarts("c1ccccc1O[C](=O)")  # direct ester linkage
    ]

    for pattern in phenyl_acetate_patterns:
        if mol.HasSubstructMatch(pattern):
            return (True, "Contains phenyl acetate ester pattern")
    
    return (False, "Does not match phenyl acetate ester pattern")