"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid has a specific steroid backbone (four fused rings) with a ketone group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for the steroid backbone - four fused rings
    # This includes the appropriate ring fusions: cyclohexane-cyclohexane-cyclohexane-cyclopentane
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C2C1)CCC4C3CC(C4)C")
    
    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"

    # SMARTS pattern for a 3-oxo group in the steroid backbone, considering it's on the A-ring
    oxo_group_pattern = Chem.MolFromSmarts("C1([#6])-C(=O)-[!#1]-[!#1]-1")  # Pattern refinement might be required

    # Identify the presence of the oxo group at position 3 on a steroid backbone
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo group found at position 3"
    
    return True, "Contains steroid backbone with 3-oxo group"

# Note: Additional structure-specific refinements may be necessary based on further testing and examples.