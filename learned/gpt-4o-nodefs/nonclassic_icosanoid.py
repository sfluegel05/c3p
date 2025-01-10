"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    Nonclassic icosanoids are characterized by a mix of polyunsaturated chains,
    potential epoxy and hydroxyl groups, and typically a carboxylic acid group.

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
    
    # Generalized pattern for any epoxy form (including larger rings)
    epoxy_pattern = Chem.MolFromSmarts("C-O-C")
    if not mol.HasSubstructMatch(epoxy_pattern):
        return False, "No generalized epoxy group found"
    
    # Looser polyunsaturated chain (at least three double bonds, not consecutive)
    poly_pattern = Chem.MolFromSmarts("C=C.[C=C]{2,}")
    if not mol.HasSubstructMatch(poly_pattern):
        return False, "No loosely polyunsaturated chain pattern found"

    # Count hydroxyl groups (loosening stereochemistry)
    hydroxyl_pattern = Chem.MolFromSmarts("C(O)")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl group(s), require at least 2"

    # Check for carboxylic acid group presence
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Contains generalized epoxy group, loose polyunsaturated chain, sufficient hydroxyl groups, and carboxylic acid"