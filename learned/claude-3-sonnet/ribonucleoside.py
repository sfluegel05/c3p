"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside consists of a nucleobase attached to a D-ribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ribose with correct stereochemistry
    # Beta-D-ribofuranose pattern allowing for modifications
    ribose_pattern = Chem.MolFromSmarts("[CH2X4][C@H]1O[C@H]([*])[C@H]([OH1,OH0])[C@@H]1[OH1,OH0]")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No beta-D-ribose sugar found"

    # Basic nucleobase patterns (more permissive)
    base_patterns = [
        # Core patterns
        "[n]1[c]([*,#1,N])[n,c][c,n][c,n]1",  # Basic azine ring
        "[n]1[c]([*,#1,N])[c]([*,#1,N])[c,n][c,n]1", # Substituted azine
        "[n]1[c][n,c][c]([*,#1,N])[c,n]1",  # Another azine variant
        
        # Pyrimidines and variants
        "[nX3]1[c]([*,#1,=O,=S,N])[n,c][c,n][c,n]1",  # Generic pyrimidine
        "[nX3]1[c]([*,#1,=O,=S,N])[n][c]([*,#1,=O,=S,N])[c,n]1", # Modified pyrimidine
        
        # Purines and variants
        "[c,n]1[n][c]2[c]([n,c]1)[n,c][c,n][n,c]2",  # Basic purine
        "[c,n]1[n][c]2[c]([n,c]1)[n,c]([*,#1,N])[c,n][n,c]2", # Modified purine
        
        # Additional patterns for modified bases
        "[n]1[c]([C,N,O,S])[n,c][c,n][c,n]1",  # Substituted base
        "[n]1[c]2[n][c]([*,#1,N])[n,c][c]2[n,c]1", # Fused system
    ]

    has_base = False
    for pattern in base_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_base = True
            break
            
    if not has_base:
        return False, "No recognized nucleobase found"

    # Check for N-glycosidic bond (more permissive)
    n_glycosidic_pattern = Chem.MolFromSmarts("[#7][C@H]1O[C@H]([CH2][OH1,OH0])[C@H]([OH1,OH0])[C@@H]1[OH1,OH0]")
    if not mol.HasSubstructMatch(n_glycosidic_pattern):
        return False, "No proper N-glycosidic bond found"

    # Ring count check
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count"

    # Success case
    return True, "Contains beta-D-ribose sugar connected to nucleobase via N-glycosidic bond"