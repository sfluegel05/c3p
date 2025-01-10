"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:28874 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    A 1-phosphatidyl-1D-myo-inositol consists of a glycerol backbone with two fatty acid chains,
    a phosphate group, and a 1D-myo-inositol moiety attached to the phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the complete pattern for 1-phosphatidyl-1D-myo-inositol
    pattern = Chem.MolFromSmarts(
        "[CH2X4][CHX4]([OX2][CX3](=[OX1]))[CH2X4][OX2]P(=O)([OX2])[OX2][C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@@H]1O)O)O)O)O)"
    )
    
    # Check if the complete pattern matches
    if not mol.HasSubstructMatch(pattern):
        return False, "Does not match 1-phosphatidyl-1D-myo-inositol pattern"

    # Verify the presence of two fatty acid chains
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2][CX3](=[OX1])"))
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Verify the stereochemistry of the inositol moiety
    inositol_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@@H]1O)O)O)O)O)")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "Incorrect stereochemistry in inositol moiety"

    # Verify the connectivity between components
    connectivity_pattern = Chem.MolFromSmarts(
        "[CH2X4][CHX4]([OX2][CX3](=[OX1]))[CH2X4][OX2]P(=O)([OX2])[OX2][C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@@H]1O)O)O)O)O)"
    )
    if not mol.HasSubstructMatch(connectivity_pattern):
        return False, "Incorrect connectivity between components"

    # Verify the phosphate group is properly positioned
    phosphate_pattern = Chem.MolFromSmarts("[CH2X4][OX2]P(=O)([OX2])[OX2][C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@@H]1O)O)O)O)O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not properly positioned"

    # Additional check for the correct number of hydroxyl groups on the inositol moiety
    hydroxyl_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@@H]1O)O)O)O)O)"))
    if len(hydroxyl_matches) != 5:
        return False, "Incorrect number of hydroxyl groups on the inositol moiety"

    return True, "Contains glycerol backbone with two fatty acid chains, a phosphate group, and a 1D-myo-inositol moiety attached to the phosphate group"