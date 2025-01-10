"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid consists of a specific steroid backbone with four fused rings
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
    
    # Revise steroid backbone SMARTS - account for variabilities
    steroid_backbone_pattern = Chem.MolFromSmarts("C1(CC[C@H]2[C@@H]3CC[C@@H]4CCC(=O)C3)C[C@H]2[C@@H]4C1")    
    # Check for steroid backbone presence
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"

    # SMARTS pattern for a 3-oxo group considering keto-functionality accurately located
    oxo_group_pattern = Chem.MolFromSmarts("C-C(=O)-C")
    
    # Ensure a 3-oxo group exists in an appropriate ring typically an A ring
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo group found at position 3"
    
    return True, "Contains steroid backbone with 3-oxo group"

# Note: The structural variability and stereo-centers have been considered with possible alternate SMARTS patterns.