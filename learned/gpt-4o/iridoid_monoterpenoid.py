"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    The molecule should contain a cyclopentane ring fused to a six-membered oxygen heterocycle,
    with variations observed in the class of iridoid monoterpenoids.

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

    # Check patterns indicative of iridoid monoterpenoids
    # Expanded pattern to capture stereochemistry and common ring configurations
    iridoid_pattern = Chem.MolFromSmarts("C1[C@H]CC([C@@H]2O[C@@H]2)=O")  # Adapted from common iridoid core structures
    secoiridoid_pattern = Chem.MolFromSmarts("O1C=CC(O)[C@H]1")  # Representing common secoiridoid fragment structures

    is_iridoid = mol.HasSubstructMatch(iridoid_pattern)
    is_secoiridoid = mol.HasSubstructMatch(secoiridoid_pattern)

    if is_iridoid:
        return True, "Matches an iridoid structure with cyclopentane fused to a six-membered oxygen heterocycle."
    
    if is_secoiridoid:
        return True, "Matches a secoiridoid structure, characterized by specific ring opening."

    return False, "Does not match characteristic iridoid monoterpenoid structures"