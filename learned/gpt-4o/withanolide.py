"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    Withanolides are C28 steroid lactones with modified side chains forming lactone rings and substituted derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Check for the steroid backbone (a more generic pattern for four fused rings)
    steroid_pattern = Chem.MolFromSmarts("C1C[C@@H]2[C@@H](CC[C@]3([C@@]2(CC[C@H]4[C@@]3(CCCC4)C)C)C)C1")  # Upgraded SMARTS for steroids
    if not mol.HasSubstructMatch(steroid_pattern):
        return (False, "No steroid backbone detected")
    
    # Flexible pattern for lactone rings
    lactone_pattern = Chem.MolFromSmarts("C1OC(=O)C=CC1|C1O[C@H](C)C(=O)C=C1")  # Consider different lactone ring sizes/patterns
    if not mol.HasSubstructMatch(lactone_pattern):
        return (False, "No lactone group found")
    
    # Broad feature detection for functional groups (hydroxyl, ketone)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    ketone_pattern = Chem.MolFromSmarts("C(=O)[C]")
    if not (mol.HasSubstructMatch(hydroxyl_pattern) or mol.HasSubstructMatch(ketone_pattern)):
        return (False, "No hydroxyl or ketone functionalities detected")

    # Success case
    return (True, "Contains features consistent with withanolides: steroid backbone, lactone ring, and functional groups")