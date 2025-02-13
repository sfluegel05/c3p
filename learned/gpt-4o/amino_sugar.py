"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is defined as any sugar having one or more alcoholic hydroxy groups
    replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Recognize sugar-like rings (adjusted to catch broader structures)
    sugar_patterns = [
        Chem.MolFromSmarts("[O;R1]1[C;R1][C;R1][C;R1][C;R1][C;R1][C;R1]1"),  # Generic 6-membered ring
        Chem.MolFromSmarts("[O;R1]1[C;R1][C;R1][C;R1][C;R1][O;R1]1"),        # 5-membered ring (furanose-like)
    ]
    
    found_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(pattern):
            found_sugar = True
            break
    
    if not found_sugar:
        return False, "No typical sugar-like ring structure found"

    # Detect amino group replacing OH
    amino_group_pattern = Chem.MolFromSmarts("[NX3H2,NX3H,NX3](C)-[OH]")
    if not mol.HasSubstructMatch(amino_group_pattern):
        return False, "No amino group replacing a hydroxy group found"
    
    return True, "Contains a sugar-like ring structure with one or more hydroxyl groups replaced by amino groups"