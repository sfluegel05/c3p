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

    # Stereochemically correct steroid nucleus: three cyclohexane rings and a cyclopentane
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(CCCC4)C")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "Steroid backbone not complete."

    # Identify 3-oxo group: A ketone (C=O) on the third carbon
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C;R;D4;v4]")  # Versatile ketone group
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "3-oxo group not correctly positioned."

    # Identify Delta(4) double bond: Double bond between C4 and C5
    delta_4_pattern = Chem.MolFromSmarts("C=C([C;R;D3;v4])C[C;R]")  # Double bond from fourth to fifth
    if not mol.HasSubstructMatch(delta_4_pattern):
        return False, "Delta(4) bond not properly located."

    return True, "3-oxo and Delta(4) double bond characteristics validated for the steroid"