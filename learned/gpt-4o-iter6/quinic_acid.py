"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or its derivative based on its SMILES string.
    
    Quinic acid is a cyclitol carboxylic acid derivative characterized by a cyclohexane ring
    with multiple hydroxyl groups and a carboxylic acid group.
    The molecule may have additional ester-linked groups like caffeoyl or feruloyl.

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

    # More generalized pattern for quinic acid backbone
    quinic_acid_pattern = Chem.MolFromSmarts("C1(CC(C(C(C1O)O)O)O)C(=O)O")
    
    # Look for stereochemistry-independent quinic acid pattern
    if not mol.HasSubstructMatch(quinic_acid_pattern):
        return False, "Missing quinic acid backbone structure"

    # Patterns for common ester-linked modifications
    caffeoyl_pattern = Chem.MolFromSmarts("c1ccc(O)c(O)c1/C=C/C(=O)O")
    feruloyl_pattern = Chem.MolFromSmarts("c1cc(O)c(OC)cc1/C=C/C(=O)O")
    
    # General ester pattern linked to cyclohexane ring
    general_ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]1C(O)C(O)C(O)C(C1)O")

    caffeoyl_matches = mol.HasSubstructMatch(caffeoyl_pattern)
    feruloyl_matches = mol.HasSubstructMatch(feruloyl_pattern)
    ester_matches = mol.HasSubstructMatch(general_ester_pattern)

    if caffeoyl_matches:
        return True, "Quinic acid derivative with caffeoyl group(s)"
    if feruloyl_matches:
        return True, "Quinic acid derivative with feruloyl group(s)"
    if ester_matches:
        return True, "Quinic acid derivative with esterified oxy groups"

    # Default to basic quinic acid if only backbone is present
    return True, "Base quinic acid present"