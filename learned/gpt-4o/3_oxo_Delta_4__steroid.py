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
        bool: True if the molecule is a 3-oxo-Δ(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broader and flexible steroid backbone pattern
    # Allows for steroid core: [C@@]1CC[C@]2[C@@]3C=CC(=O)CC[C@]3(C)CC[C@@]2(C)[C@H]1
    steroid_backbone_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4(C3)C=C[C@H]4')
    
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"

    # 3-oxo group specific to C3 position with flexibility
    oxo_group_pattern = Chem.MolFromSmarts('O=C[C;R2,C3]')

    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo group found"
    
    # Δ(4) double bond between C4 and C5
    delta_4_pattern = Chem.MolFromSmarts('C=CC(=O)[C;R3,C5]')

    if not mol.HasSubstructMatch(delta_4_pattern):
        return False, "No Δ(4) double bond found within steroid context"

    return True, "Contains steroid backbone with a 3-oxo group and a Δ(4) double bond"