"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    Nonclassic icosanoids typically feature polyunsaturated carbon chains, multiple
    hydroxyl groups, one or more epoxy groups, and a terminal carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Specific epoxy group e.g., in a 3-membered ring
    epoxy_pattern = Chem.MolFromSmarts("[C@H]1O[C@@H]1")
    if not mol.HasSubstructMatch(epoxy_pattern):
        return False, "No specific epoxy group found"
    
    # Polyunsaturated chain with at least 3 double bonds
    poly_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if not mol.HasSubstructMatch(poly_pattern):
        return False, "No suitable polyunsaturated chain found"

    # Check for sufficient hydroxyl groups, not necessarily considering stereochemistry
    hydroxyl_pattern = Chem.MolFromSmarts("CO")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl group(s), require at least 2"

    # Check for a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Contains epoxy group, polyunsaturated chain, sufficient hydroxyl groups, and carboxylic acid"