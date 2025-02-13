"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate derivative based on its SMILES string.
    A phenyl acetate can have an acetyl group (-C(=O)) connected via oxygen or nitrogen to a phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for acetate ester connected to phenyl group pattern
    # This pattern includes esters, where the carbonyl carbon is directly bonded to an oxygen attached to a phenyl ring.
    phenyl_acetate_pattern = Chem.MolFromSmarts("c1ccccc1OC(=O)C")
    
    if mol.HasSubstructMatch(phenyl_acetate_pattern):
        return True, "Contains phenyl acetate ester pattern"

    return False, "Does not match phenyl acetate ester pattern"