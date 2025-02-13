"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or its derivative based on its SMILES string.
    
    Quinic acid is a cyclitol carboxylic acid with specific hydroxyl and carboxylic acid groups on a cyclohexane ring.

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

    # Define SMARTS pattern for cyclohexane with hydroxyls and carboxylic acid groups
    quinic_acid_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)[C@H](O)[C@H](C(=O)O)C1")
    # Check for quinic acid backbone
    if not mol.HasSubstructMatch(quinic_acid_pattern):
        return False, "Missing quinic acid backbone structure"
    
    # Check for additional substitutions like caffeoyl or feruloyl groups
    caffeoyl_pattern = Chem.MolFromSmarts("c1ccc(O)c(O)c1/C=C/C(=O)O")
    feruloyl_pattern = Chem.MolFromSmarts("c1cc(O)cc(OC)c1/C=C/C(=O)O")
    
    caffeoyl_matches = mol.GetSubstructMatches(caffeoyl_pattern)
    feruloyl_matches = mol.GetSubstructMatches(feruloyl_pattern)

    if len(caffeoyl_matches) > 0:
        return True, "Quinic acid derivative with caffeoyl group(s)"
    if len(feruloyl_matches) > 0:
        return True, "Quinic acid derivative with feruloyl group(s)"

    # If no additional groups, classify as base quinic acid
    return True, "Base quinic acid with no additional groups detected"