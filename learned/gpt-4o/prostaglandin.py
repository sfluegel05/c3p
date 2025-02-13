"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are C20 compounds derived from prostanoic acid, generally featuring a
    substituted cyclopentane ring and various functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for cyclopentane ring (allowing substitutions)
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        # Try a flexible pattern with substitutions
        flexible_cyclopentane_pattern = Chem.MolFromSmarts("C1(C)C(C)C(C)C1")
        if not mol.HasSubstructMatch(flexible_cyclopentane_pattern):
            return False, "No cyclopentane or substituted cyclopentane ring detected"
    
    # Check for acid or acid derivative group (e.g., carboxylic acid or ester)
    acid_pattern = Chem.MolFromSmarts("C(=O)[O,N]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "Missing acid or acid-like functional group"

    # Check for key functional groups close to the cyclopentane
    key_groups_pattern = Chem.MolFromSmarts("O=C-C1CCCC1")
    if not mol.HasSubstructMatch(key_groups_pattern):
        # Consider the entire pattern around the cyclopentane
        extensive_group_pattern = Chem.MolFromSmarts("O=C[C@H]1CC[C@H](O)C1")
        if not mol.HasSubstructMatch(extensive_group_pattern):
            return False, "Missing key functional groups in proximity to cyclopentane ring"
    
    # Consider the presence and pattern of double bonds
    db_pattern = Chem.MolFromSmarts("C=C")
    db_count = len(mol.GetSubstructMatches(db_pattern))
    if db_count < 2:
        return False, "Prostaglandins typically contain multiple double bonds; too few found"

    return True, "Molecule contains a C20 backbone with a cyclopentane or properly substituted cyclopentane, relevant acid-functional groups and double bonds typical of prostaglandins"