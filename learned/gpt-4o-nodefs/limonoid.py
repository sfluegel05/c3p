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
    for common structural motifs seen in limonoids like highly oxygenated
    triterpene derivatives with diverse ring systems and often with a furan ring.

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

    # Flexible ring skeleton pattern (allow more variation in rings)
    flexible_core_pattern = Chem.MolFromSmarts("C1CCC2CCCC(C=C3)CC3(O)OC2C1")  # Example abstraction
    
    # Check for patterns such as furan rings which might be common but not mandatory
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")

    # Verify presence of flexible limonoid core
    if not mol.HasSubstructMatch(flexible_core_pattern):
        return False, "No flexible limonoid core structure found"

    # Check for presence of a furan ring (optional but common)
    if mol.HasSubstructMatch(furan_pattern):
        furan_presence = " with furan ring"
    else:
        furan_presence = ", furan ring not detected"

    # Count multiple oxygens: limonoids are typically highly oxygenated
    oxygen_pattern = Chem.MolFromSmarts("[O]")
    oxy_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxy_matches) < 5:
        return False, f"Found {len(oxy_matches)} oxygenated groups, expecting more than 5"

    return True, f"Contains limonoid core structure{furan_presence} and high oxygenation"