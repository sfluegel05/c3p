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

    # Pattern for 2-oxo monocarboxylic acid:
    # - Carbon with single ketone (not part of carboxyl)
    # - Connected to carbon of carboxylic acid
    # - Not part of an ester or other carboxylic acid
    alpha_keto_pattern = Chem.MolFromSmarts("[#6]-[CX3](=[OX1])-[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(alpha_keto_pattern):
        return False, "No alpha-keto acid pattern found"

    # Additional checks to exclude false positives
    
    # Exclude cases where the ketone is part of a carboxylic acid or ester
    problematic_pattern = Chem.MolFromSmarts("[OX2H1,OX2R]-[CX3](=[OX1])-[CX3](=[OX1])[OX2H1]")
    if mol.HasSubstructMatch(problematic_pattern):
        return False, "The ketone is part of another acid or ester group"

    # Check that carbon between ketone and acid is not part of other functional groups
    connecting_c_pattern = Chem.MolFromSmarts("[#6]-[CX3](=[OX1])-[CX3](=[OX1])[OX2H1]")
    matches = mol.GetSubstructMatches(connecting_c_pattern)
    
    for match in matches:
        connecting_carbon = mol.GetAtomWithIdx(match[2])
        if connecting_carbon.GetDegree() > 2:  # Should only connect to C=O and COOH
            continue
        return True, "Contains a valid 2-oxo monocarboxylic acid pattern"

    return False, "The structure does not match required 2-oxo monocarboxylic acid pattern"