"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide typically includes a modified steroid skeleton with additional
    oxygen functionalities, often in the form of a lactone group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a more flexible pattern for the steroid core
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CC(C(C4)C3)C2C1")  # Simplified 6-6-6-5 fused ring system
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid-like structure found"
    
    # Look for potential lactone groups ring closure patterns
    lactone_pattern1 = Chem.MolFromSmarts("O=C1COC=C1")  # 5-membered lactone
    lactone_pattern2 = Chem.MolFromSmarts("O=C1COCC=C1") # 6-membered lactone
    if not (mol.HasSubstructMatch(lactone_pattern1) or mol.HasSubstructMatch(lactone_pattern2)):
        return False, "No lactone ring found"

    # Check for oxygen functionalities - allow a more general check
    oxy_gen_patterns = [
        Chem.MolFromSmarts("O"),  # Oxygen atoms
        Chem.MolFromSmarts("[OX2H]"),  # Hydroxyl oxygen
        Chem.MolFromSmarts("[#6]=[OX1]"),  # Carbonyl oxygen
    ]

    oxy_count = sum(len(mol.GetSubstructMatches(patt)) for patt in oxy_gen_patterns)
    # Ensure structure has sufficient oxygen functional groups commonly seen in Withanolides
    if oxy_count < 3:
        return False, "Insufficient oxygen functionalities detected"

    return True, "Matches characteristics of withanolide: steroid core, lactone presence, and adequate oxygen functionalities"