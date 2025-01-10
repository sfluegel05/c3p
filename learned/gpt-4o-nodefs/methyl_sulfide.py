"""
Classifies: CHEBI:86315 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide contains the CH3-S motif within its structure.

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
    
    # Look for a methylthio group pattern (CH3-S)
    # This pattern considers that S atom connects to any groups including a CH3 group.
    methylthio_pattern = Chem.MolFromSmarts("[CH3]S")
    
    # Check for substructure match
    if mol.HasSubstructMatch(methylthio_pattern):
        return True, "Contains methylthio group (CH3-S)"

    return False, "No suitable methylthio group found"