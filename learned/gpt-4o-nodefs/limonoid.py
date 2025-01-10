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

    # Revise triterpene backbone pattern: Using a generic steroid-like core pattern,
    # since limonoids are derived from tetracyclic triterpenes
    triterpene_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4(C)C(=O)CCC4(C)C3CCC2C1")
    
    # Common furan pattern
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")

    # Check for broader triterpene-like backbone
    if not mol.HasSubstructMatch(triterpene_pattern):
        return False, "No triterpene or steroid-like backbone found"
    
    # Check for the presence of a furan ring
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Count oxygenation: Look for carbonyls, hydroxyls, esters, etc.
    # Allow a slightly flexible range expecting more than 5 oxygen-related functionalities for limonoids
    oxygen_pattern = Chem.MolFromSmarts("[OX2H,OX1=]")
    oxy_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxy_matches) < 5:
        return False, f"Found {len(oxy_matches)} oxygenated groups, need more than 5"

    return True, "Contains limonoid-specific structural features including a triterpene/steroid core, furan ring, and multiple oxygens"