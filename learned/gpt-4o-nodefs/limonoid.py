"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: Limonoid
"""
from rdkit import Chem

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string by checking
    for common limonoid structural motifs including tetracyclic frameworks,
    furan rings, and high oxygenation levels.

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

    # Broader tetracyclic or triterpenoid core pattern
    # Using a more abstract triterpene backbone that accommodates flexibility in ring fusion
    tetracyclic_pattern = Chem.MolFromSmarts("C1CCC2C3CCCC4CCC(CCC4)C3CCC2C1")
    
    # Common substructure featured in limonoids
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")

    # Verify tetracyclic or similar backbone
    if not mol.HasSubstructMatch(tetracyclic_pattern):
        return False, "No recognizably flexible tetracyclic backbone found"

    # Check for presence of a furan ring
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Count multiple oxygens: limonoids are highly oxygenated
    oxygen_pattern = Chem.MolFromSmarts("[O]")
    oxy_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxy_matches) < 5:
        return False, f"Found {len(oxy_matches)} oxygenated groups, expecting more than 5"

    return True, "Contains limonoid-specific tetracyclic framework, furan ring, and high oxygenation"