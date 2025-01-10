"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or its derivative based on its SMILES string.
    
    Quinic acid is a cyclitol carboxylic acid derivative characterized by a cyclohexane ring
    with multiple hydroxyl groups and a carboxylic acid group.
    The molecule may have additional ester-linked modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is quinic acid or its derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More flexible core pattern for quinic acid
    # Basic cyclohexane with optional stereochemistry, at least 3 OH and 1 COOH
    core_pattern = Chem.MolFromSmarts("C1C(C(O)C(C(C1[OH])))[OH]C(=O)O")
    
    # Check for quinic acid core
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing quinic acid backbone structure"
    
    # Pattern for caffeoyl or other ester-linked modifications
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]1C[C@H](O)C(O)C[C@H]1O")
    if mol.HasSubstructMatch(ester_pattern):
        return True, "Quinic acid derivative with ester linkages detected"
    
    # Handling additional structural complexity
    additional_complexity_pattern = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(additional_complexity_pattern):
        return True, "Quinic acid with complex modifications including acetylation"

    # If structure does not match any complex patterns but matches core
    return True, "Base quinic acid or simple derivative identified"