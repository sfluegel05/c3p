"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    Iridoid monoterpenoids typically involve a cyclopentane ring fused to a six-membered oxygen heterocycle,
    or an open lactone form as seen in secoiridoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Recognize major motif: cyclopentane ring fused with oxygen heterocycle
    cyclopentane_oxygen_pattern = Chem.MolFromSmarts("[C]1[C][C][C][O]1")
    if mol.HasSubstructMatch(cyclopentane_oxygen_pattern):
        # Check for variant patterns indicating iridoid structures
        cyclopentane_fused_pattern = Chem.MolFromSmarts("[C]1[C][C][O][C]2[C][C]1O2")
        if mol.HasSubstructMatch(cyclopentane_fused_pattern):
            return True, "Found cyclopentane fused with six-membered oxygen-containing heterocycle"

    # Recognize open lactone structures indicative of secoiridoids
    secoiridoid_pattern = Chem.MolFromSmarts("[C]1(=O)[O][C][C]=C[C][O1]")
    if mol.HasSubstructMatch(secoiridoid_pattern):
        return True, "Open lactone structure typical of secoiridoids detected"

    # Check for additional patterns of iridoid variations (e.g., complex lactones)
    iridoid_complex_lactone_pattern = Chem.MolFromSmarts("[C]1[C]([C][C])C[C]OC1")
    if mol.HasSubstructMatch(iridoid_complex_lactone_pattern):
        return True, "Complex iridoid structural motif with lactone content recognized"

    return False, "Iridoid monoterpenoid characteristic structure not detected"