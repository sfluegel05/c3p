"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: CHEBI:75768 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    Neoflavonoids have a 1-benzopyran core with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1-benzopyran core pattern with position 4 marked
    # The [#0] is a dummy atom marking position 4
    benzopyran_pattern = Chem.MolFromSmarts("O1c2ccccc2C([#0])C=C1")
    
    # Check if the basic core is present
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 1-benzopyran core found"

    # Pattern for 1-benzopyran with aryl group at position 4
    # More specific pattern that requires the aryl group to be directly attached
    neoflavonoid_pattern = Chem.MolFromSmarts("O1c2ccccc2C(c3ccccc3)C=C1")
    
    if not mol.HasSubstructMatch(neoflavonoid_pattern):
        return False, "No aryl substituent at position 4"

    # Additional check for variations in the aryl group
    # This allows for substituted phenyl rings
    extended_aryl_pattern = Chem.MolFromSmarts("O1c2ccccc2C(c3[c,n]cccc3)C=C1")
    
    # Count matches to determine confidence
    basic_matches = len(mol.GetSubstructMatches(neoflavonoid_pattern))
    extended_matches = len(mol.GetSubstructMatches(extended_aryl_pattern))
    
    confidence = "medium"
    if basic_matches > 0 and extended_matches > 0:
        confidence = "high"

    # Common substituents in neoflavonoids
    substituents = [
        (Chem.MolFromSmarts("[OH]"), "hydroxyl"),
        (Chem.MolFromSmarts("O[CH3]"), "methoxy"),
        (Chem.MolFromSmarts("C(=O)"), "carbonyl")
    ]
    
    found_substituents = []
    for pattern, name in substituents:
        if mol.HasSubstructMatch(pattern):
            found_substituents.append(name)
    
    # Build detailed reason
    reason = "Contains 1-benzopyran core with aryl substituent at position 4"
    if found_substituents:
        reason += f" and {', '.join(found_substituents)} groups"
    reason += f" ({confidence} confidence)"

    return True, reason