"""
Classifies: CHEBI:86315 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    Methyl sulfides contain a methyl group attached to a sulfur atom (CH3-S-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for methylthio group pattern (C-S)
    methylthio_pattern = Chem.MolFromSmarts("CS")
    if mol.HasSubstructMatch(methylthio_pattern):
        return True, "Contains methylthio group (CH3-S-)"
    
    return False, "No methylthio group found"