"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: CHEBI:37783 3-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid has a steroid backbone with a ketone group at position 3.

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

    # Define a more flexible steroid backbone pattern (four fused rings)
    # This pattern captures the four-fused-ring structure with possible double bonds and substitutions
    steroid_backbone_pattern = Chem.MolFromSmarts("[C@]12[C@]3([C@]4([C@]([C@H]1CC2)(CC3)CC4)C)")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        # Try a more general pattern if the first one fails
        steroid_backbone_pattern = Chem.MolFromSmarts("[C@]12[C@]3([C@]4([C@]([C@H]1CC2)(CC3)CC4)C)")
        if not mol.HasSubstructMatch(steroid_backbone_pattern):
            return False, "No steroid backbone found"

    # Define the 3-oxo pattern (ketone at position 3)
    # The ketone should be attached to the carbon at position 3 in the steroid backbone
    oxo_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[C@]12[C@]3([C@]4([C@]([C@H]1CC2)(CC3)CC4)C)")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No ketone group found at position 3"

    return True, "Contains a steroid backbone with a ketone group at position 3"