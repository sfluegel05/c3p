"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: CHEBI:17971 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def has_ketone_at_position_3(mol):
    """
    Checks if the molecule has a ketone group at position 3.
    """
    ketone_pattern = Chem.MolFromSmarts("[C@H]1[C@H](CC[C@@]2([C@]1([C@H]([C@@]3([C@H](CC2)C)C)(C)C)C)C(=O)")
    if ketone_pattern is None:
        return False
    return mol.HasSubstructMatch(ketone_pattern)

def has_beta_config_at_position_5(mol):
    """
    Checks if the molecule has a beta configuration at position 5.
    """
    beta_pattern_1 = Chem.MolFromSmarts("[C@@]1([C@H](CC[C@@]2([C@]1([C@H]([C@@]3([C@H](CC2)C)C)(C)C)C)C)C")
    beta_pattern_2 = Chem.MolFromSmarts("[C@@]12([C@H](CC[C@@]3([C@]1([C@H]([C@@]4([C@H](CC3)C)C)(C)C)C)C)CC[C@@H]2C")
    if beta_pattern_1 is None or beta_pattern_2 is None:
        return False
    return mol.HasSubstructMatch(beta_pattern_1) or mol.HasSubstructMatch(beta_pattern_2)

def is_steroid_backbone(mol):
    """
    Checks if the molecule has a steroid backbone.
    """
    steroid_pattern_1 = Chem.MolFromSmarts("[C@]12CC[C@@]3([C@@]1(CCC[C@@H]2O)C)[C@H](CC[C@@H]4[C@]3(CCC(=O)C[C@@H]4)C)C")
    steroid_pattern_2 = Chem.MolFromSmarts("[C@]12CC[C@@]3([C@@]1(CCC[C@@H]2O)C)[C@H](CC[C@@H]4[C@]3(CC[C@@H](C[C@@H]4)O)C)C")
    if steroid_pattern_1 is None or steroid_pattern_2 is None:
        return False
    return mol.HasSubstructMatch(steroid_pattern_1) or mol.HasSubstructMatch(steroid_pattern_2)

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is a steroid with a ketone at position 3 and beta configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule has all the required features
    has_ketone_3 = has_ketone_at_position_3(mol)
    has_beta_5 = has_beta_config_at_position_5(mol)
    is_steroid = is_steroid_backbone(mol)

    if has_ketone_3 and has_beta_5 and is_steroid:
        return True, "Molecule contains a ketone at position 3 and beta configuration at position 5 on a steroid backbone"
    elif not has_ketone_3:
        return False, "No ketone group found at position 3"
    elif not has_beta_5:
        return False, "No beta configuration found at position 5"
    else:
        return False, "Not a steroid backbone"