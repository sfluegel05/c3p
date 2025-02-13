"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are a type of flavonoid characterized by a 2-phenyl-3,4-dihydro-2H-chromen-3-ol structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Catechin pattern: 2-phenyl-3,4-dihydro-2H-chromen-3-ol
    catechin_pattern = Chem.MolFromSmarts('c1cc(c2c(c1)CC(O)c3c2c(ccc3)O)O')
    if not mol.HasSubstructMatch(catechin_pattern):
        return False, "Does not match catechin core structure (2-phenyl-3,4-dihydro-2H-chromen-3-ol)"
    
    # Check for multiple hydroxyl groups on aromatic rings that is characteristic of catechins
    hydroxyl_pattern = Chem.MolFromSmarts('c[cH][cH][O]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    min_hydroxyl_count = 4  # Catechins commonly have multiple hydroxyl groups
    if len(hydroxyl_matches) < min_hydroxyl_count:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need at least {min_hydroxyl_count}"

    return True, "Contains catechin core structure with sufficient hydroxylation"

# Example of how you might use this function:
smiles_example = "O1[C@@H]([C@@H](O)CC2=C1C=C(O)C=C2)C3=CC(O)=C(O)C=C3"  # (-)-catechin
result, reason = is_catechin(smiles_example)
print(f"Result: {result}, Reason: {reason}")