"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid consists of a steroid backbone characterized by four fused rings
    and a ketone group at position 3 on the A-ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Revise steroid backbone SMARTS - account for four fused rings without specifying stereochemistry detailedly
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2CCC3C4CCCC(C4)C3C2C1")    
    
    # Check for steroid backbone presence
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"

    # Improved SMARTS pattern for a 3-oxo group on the A ring of the steroid
    oxo_group_pattern = Chem.MolFromSmarts("C1(=O)CC2C3CCCC3CC12")  # Ketone at ring joining
    
    # Ensure a 3-oxo group exists in an appropriate location (A-ring)
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo group found at position 3"
    
    return True, "Contains steroid backbone with 3-oxo group"

# Note: Adjustments have been made to better reflect 3-oxo-ketone position within the A-ring and improve steroid backbone recognition.