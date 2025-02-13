"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid is characterized by a 3-oxo group and a double bond 
    between carbon 4 and 5 of the steroid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid backbone: 3 cyclohexane rings and 1 cyclopentane ring fused together
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CC2CC3CC4CCC(C3)C2C=C1C4")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "Steroid backbone not complete."

    # 3-Oxo group: a keto group attached at the C3 position
    oxo_pattern = Chem.MolFromSmarts("C1(C=O)=CC(C)([C;R2])C=C1")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "3-oxo group not correctly positioned."

    # Delta(4) double bond: C=C bond between the 4th and 5th carbon of the steroid nucleus
    delta_4_pattern = Chem.MolFromSmarts("C=C(C)C1=C")
    if not mol.HasSubstructMatch(delta_4_pattern):
        return False, "Delta(4) bond not properly located."

    return True, "3-oxo and Delta(4) double bond characteristics validated for the steroid"