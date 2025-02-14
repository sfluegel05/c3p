"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: CHEBI:25495 phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid is a derivative of glycerol where one hydroxyl group is esterified with phosphoric acid and the other two with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is phosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more general glycerol backbone pattern with flexible ester groups
    glycerol_pattern = Chem.MolFromSmarts("C(CO)O")
    
    # Define phosphoric acid ester pattern, allow flexibility
    phosphoester_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    
    # Phosphatidic acid must match the flexible glycerol with esterified positions
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No flexible glycerol backbone found"
    
    # Ensure presence of the phosphoric acid ester group
    if not mol.HasSubstructMatch(phosphoester_pattern):
        return False, "No phosphoric acid ester group found"

    # Check for two esterified fatty acid groups, allowing flexibility
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    # specifically check for the length of fatty acids to be more than 10 carbons
    fatty_acid_count = sum(1 for match in ester_matches if len(match) > 10)
    if fatty_acid_count < 2:
        return False, "Less than 2 fatty acid ester groups found"

    return True, "Contains glycerol backbone with 2 fatty acid esters and a phosphoric acid ester"