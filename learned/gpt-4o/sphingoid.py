"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid generally has a long carbon chain, amino group, and hydroxyl groups,
    including potential double bonds and specific stereochemistry.
    
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
    
    # Check for long aliphatic chain pattern (14-30+ carbons)
    # Consider various chain lengths and potential unsaturations
    long_chain_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCCCCCC"),  # At least 14 contiguous carbons
        Chem.MolFromSmarts("CCCCCCCCCCCCCCC"), # At least 15 contiguous carbons
    ]
    if not any(mol.HasSubstructMatch(chain) for chain in long_chain_patterns):
        return False, "No sufficient long carbon chain recognized"
    
    # Look for the sphingoid core (2-amino-1,3-diol configuration)
    core_patterns = [
        Chem.MolFromSmarts("[C@@H](O)[C@@H](N)CO"),  # Common stereochemistry
        Chem.MolFromSmarts("[C@H](O)[C@H](N)CO"),
        Chem.MolFromSmarts("[C@H]([NH3+])[C@H](O)CO"),  # Protonated amino group
    ]
    if not any(mol.HasSubstructMatch(core) for core in core_patterns):
        return False, "No recognizable sphingoid core substructure"

    # Count hydroxyl groups (at least 2 in sphingoids)
    hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
    if hydroxyl_count < 2:
        return False, f"Insufficient hydroxyl groups; found {hydroxyl_count}"

    # Detect double bonds if present
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(double_bond_pattern):
        return True, "Contains typical sphingoid structure with additional features like double bonds"

    return True, "Contains typical structure of a sphingoid (long chain, amino, hydroxyl groups)"