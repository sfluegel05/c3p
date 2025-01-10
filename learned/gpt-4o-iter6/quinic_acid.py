"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or its derivative based on its SMILES string.
    
    Quinic acid is a cyclitol carboxylic acid with specific hydroxyl and carboxylic acid groups on a cyclohexane ring.
    The molecule may have additional groups like caffeoyl or feruloyl.

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

    # Define a more generalized SMARTS pattern for cyclohexane with hydroxyl and carboxylic acid groups
    quinic_acid_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@@H](O)C(O)C(O)[C@@H](O)C1C(=O)O")
    if not mol.HasSubstructMatch(quinic_acid_pattern):
        return False, "Missing quinic acid backbone structure"

    # Look for modification patterns like caffeoyl or feruloyl
    caffeoyl_pattern = Chem.MolFromSmarts("c1ccc(O)c(O)c1/C=C/C(=O)O")
    feruloyl_pattern = Chem.MolFromSmarts("c1cc(O)c(OC)cc1/C=C/C(=O)O")

    caffeoyl_matches = mol.GetSubstructMatches(caffeoyl_pattern)
    feruloyl_matches = mol.GetSubstructMatches(feruloyl_pattern)

    if len(caffeoyl_matches) > 0:
        return True, "Quinic acid derivative with caffeoyl group(s)"
    if len(feruloyl_matches) > 0:
        return True, "Quinic acid derivative with feruloyl group(s)"

    # Check for other common ester linkages using general patterns
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]1C([C@@H](O)C(=O))C(O)C[C@H]1O")
    if mol.HasSubstructMatch(ester_pattern):
        return True, "Quinic acid derivative with esterified oxy groups"

    # If nothing additional is matched, assume it's a form of base quinic acid
    return True, "Quinic acid with no recognized additional groups"