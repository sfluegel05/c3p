"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids generally have a long carbon chain, an amino group, and hydroxyl groups,
    as well as potential additional features like double bonds or certain stereochemistry.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a long aliphatic carbon chain (14-30+ carbons)
    long_chain_pattern = Chem.MolFromSmarts("C" * 14)  # Simplified pattern for long aliphatic chain
    if len(mol.GetSubstructMatches(long_chain_pattern)) < 1:
        return False, "No sufficient long carbon chain recognized"

    # Check for the presence of 2-amino-1,3-diol configuration, allowing variations
    core_pattern_variants = [
        Chem.MolFromSmarts("[C@H](O)[C@H](N)CO"),  # Original pattern
        Chem.MolFromSmarts("[C@H](N)[C@H](O)CO"),  # Variation
        Chem.MolFromSmarts("[C@H](O)[C@H](N)COC"), # For different hydroxyl connections
        Chem.MolFromSmarts("[C@H](N)C(O)CO"),  # Possible stereochemistry variation
    ]
    
    if not any(mol.HasSubstructMatch(cp) for cp in core_pattern_variants):
        return False, "No recognizable sphingoid core substructure"

    # Verify if there are appropriate hydroxyl groups numbers (i.e., at least 2)
    hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
    if hydroxyl_count < 2:
        return False, "Less than the required number of hydroxyl groups found"

    # Check for common side features like double bonds or additional functional groups
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(double_bond_pattern):
        return True, "Contains additional features like double bonds, consistent with sphingoid derivatives"

    return True, "Contains typical structure of a sphingoid (long chain, amino, hydroxyl groups)"