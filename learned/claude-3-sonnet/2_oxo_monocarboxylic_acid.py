"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: CHEBI:35755 2-oxo monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid has a single carboxylic acid group with a ketone group
    at the alpha position (2-oxo).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    
    if len(carboxylic_matches) == 0:
        return False, "No carboxylic acid group found"
    elif len(carboxylic_matches) > 1:
        return False, "Multiple carboxylic acid groups found"

    # Pattern for alpha-keto acid (2-oxo monocarboxylic acid)
    # Matches both direct connection and cases with substituents
    alpha_keto_patterns = [
        # Direct alpha-keto acid pattern
        Chem.MolFromSmarts("[CX3](=[OX1])[CX3](=[OX1])[OX2H1]"),
        # Pattern with possible substituents
        Chem.MolFromSmarts("[#6][CX3](=[OX1])[CX3](=[OX1])[OX2H1]"),
        # Alternative pattern for resonance forms
        Chem.MolFromSmarts("[#6]-[CX3](=O)-[CX3](=O)[OX2H1]")
    ]
    
    found_pattern = False
    for pattern in alpha_keto_patterns:
        if mol.HasSubstructMatch(pattern):
            found_pattern = True
            break
            
    if not found_pattern:
        return False, "No alpha-keto acid pattern found"

    # Exclude problematic cases
    
    # Exclude compounds with multiple ketone groups adjacent to the acid
    multiple_ketone_pattern = Chem.MolFromSmarts("[OX1]=[CX3]-[CX3](=[OX1])[CX3](=[OX1])[OX2H1]")
    if mol.HasSubstructMatch(multiple_ketone_pattern):
        return False, "Multiple ketone groups found adjacent to acid"

    # Exclude cases where the ketone is part of an anhydride or similar structure
    anhydride_pattern = Chem.MolFromSmarts("[OX2]-[CX3](=[OX1])[CX3](=[OX1])[OX2H1]")
    if mol.HasSubstructMatch(anhydride_pattern):
        return False, "Ketone is part of an anhydride or similar structure"

    # Check that the ketone carbon is not part of a ring system containing the acid
    ring_pattern = Chem.MolFromSmarts("[CX3]1(=[OX1])[CX3](=[OX1])[OX2H1]1")
    if mol.HasSubstructMatch(ring_pattern):
        return False, "Ketone and acid are part of the same ring"

    return True, "Contains a valid 2-oxo monocarboxylic acid pattern"