"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    Nonclassic icosanoids are derivatives characterized by polyunsaturated chains, 
    epoxy groups, and multiple hydroxyl groups.

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
    
    # Look for epoxy group (epoxide, cyclic ether: -O-)
    epoxy_pattern = Chem.MolFromSmarts("C1OC1")
    if not mol.HasSubstructMatch(epoxy_pattern):
        return False, "No epoxy group found"
    
    # Look for polyunsaturated chain (conjugated double bonds pattern)
    polyunsaturated_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if not mol.HasSubstructMatch(polyunsaturated_pattern):
        return False, "No polyunsaturated chain pattern found"

    # Look for multiple hydroxyl groups, accounting for stereochemistry
    hydroxyl_pattern = Chem.MolFromSmarts("[C@H](O)")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl group(s), require at least 2 with stereochemistry"

    # Check for carboxylic acid group presence
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Contains epoxy group, polyunsaturated chain, multiple hydroxyl groups, and carboxylic acid"