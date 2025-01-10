"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: 3-oxo-Δ(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Δ(4) steroid based on its SMILES string.
    A 3-oxo-Δ(4) steroid is defined by: A 3-oxo steroid conjugated to a C=C double bond
    at the alpha,beta position (typically between C4 and C5 on the steroid backbone).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Δ(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify 4-ring steroid backbone (ABCD ring system), a more robust pattern
    steroid_backbone_pattern = Chem.MolFromSmarts('C1CC2CCC3C4CCC(=O)C=CC4CCC3C2C1')
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"

    # Identify 3-oxo group (ketone at the C3 position)
    oxo_group_pattern = Chem.MolFromSmarts('C(=O)C[C;R]')
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo group found"
    
    # Identify Δ(4) double bond (between C4 and C5)
    delta_4_pattern = Chem.MolFromSmarts('C=C[C;R]')
    if not mol.HasSubstructMatch(delta_4_pattern):
        return False, "No Δ(4) double bond found"

    return True, "Contains steroid backbone with 3-oxo group and Δ(4) double bond"