"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:27555 flavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanone(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string using structural rules.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the core flavanone skeleton with possible substitutions
    flavanone_pattern = Chem.MolFromSmarts("[C&ring]1(=O)C=C(O[C&ring]1[C&ring]2ccccc2)[C&ring]")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Does not contain the flavanone core skeleton"

    # Check for common substituents (hydroxyl, methoxy, prenyl, etc.)
    substituents_pattern = Chem.MolFromSmarts("[OH,OC,C=CC(C)=C]")
    if not mol.HasSubstructMatch(substituents_pattern):
        return False, "No common flavanone substituents found"

    # Check for additional ring systems (fused, bridged, etc.)
    ring_systems_pattern = Chem.MolFromSmarts("[C&ring]1[C&ring][C&ring][C&ring][C&ring][C&ring]1")
    if mol.HasSubstructMatch(ring_systems_pattern):
        return True, "Contains the flavanone core with additional ring systems"

    # Check for specific stereochemistry at C2 (if present)
    chiral_center_pattern = Chem.MolFromSmarts("[C@H]1[C@@](O)(CC(=O)C=C(O1)c2ccccc2)[C&ring]")
    if mol.HasSubstructMatch(chiral_center_pattern):
        return True, "Contains the flavanone core with correct stereochemistry at C2"

    return True, "Matches the structural rules for flavanones"