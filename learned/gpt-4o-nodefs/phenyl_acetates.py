"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate derivative based on its SMILES string.
    A phenyl acetate has an acetyl group (-C(=O)C) directly attached to a phenyl ring.

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
    
    # Phenyl ring with acetyl group pattern
    # Specifically look for a bond from the acetyl carbonyl (C=O) to an aromatic carbon
    phenyl_acetate_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)C")
    
    if mol.HasSubstructMatch(phenyl_acetate_pattern):
        return True, "Contains phenyl ring with acetyl group attached"

    return False, "No acetyl group directly attached to phenyl ring"