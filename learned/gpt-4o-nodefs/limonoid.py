"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: Limonoid
"""
from rdkit import Chem

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Common limonoid pattern: Triterpene backbone with multiple oxygenated groups and a furan ring
    # Note: This is a simplified representation and may not be exhaustive
    triterpene_pattern = Chem.MolFromSmarts("C1(C)[C@]234[C@H](C[C@@]1(C)CC3=CC(=O)C4)C2")
    furan_pattern = Chem.MolFromSmarts("c1coc[cH]c1")
    oxygen_pattern = Chem.MolFromSmarts("[OX2,OX1=]")

    # Check for triterpene backbone
    if not mol.HasSubstructMatch(triterpene_pattern):
        return False, "No triterpene backbone found"
    
    # Check for presence of a furan ring
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Check for multiple oxygenated groups
    oxy_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxy_matches) < 3:
        return False, f"Found {len(oxy_matches)} oxygenated groups, need at least 3"

    return True, "Contains limonoid-specific structural features including a triterpene core, furan ring, and multiple oxygens"