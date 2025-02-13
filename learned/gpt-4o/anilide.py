"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is defined as any aromatic amide obtained by acylation of aniline.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aniline-like structure pattern (phenyl group bonded to nitrogen)
    aniline_pattern = Chem.MolFromSmarts("c1ccccc1NC(=O)")
    if mol.HasSubstructMatch(aniline_pattern):
        return True, "Contains aniline-like aromatic amide structure"
    else:
        return False, "Missing aniline-like structure characteristic of anilides"

    return False, "Does not match anilide structural criteria"