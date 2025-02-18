"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: CHEBI:xxxxx oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is a fatty acid containing at least one aldehydic or ketonic group in addition to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_matches:
        return False, "No carboxylic acid group found"

    # Collect carbonyl carbons from carboxylic acids
    carboxylic_carbonyl_carbons = {match[0] for match in carboxylic_matches}

    # Check for aldehydic groups (CH=O) not part of carboxylic acid
    aldehyde_pattern = Chem.MolFromSmarts("[CH1]=O")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    for match in aldehyde_matches:
        carbonyl_carbon = match[0]
        if carbonyl_carbon not in carboxylic_carbonyl_carbons:
            return True, "Contains aldehyde group in addition to carboxylic acid"

    # Check for ketonic groups (C-C(=O)-C) not part of carboxylic acid
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)([#6])[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    for match in ketone_matches:
        carbonyl_carbon = match[0]
        if carbonyl_carbon not in carboxylic_carbonyl_carbons:
            return True, "Contains ketone group in addition to carboxylic acid"

    return False, "No aldehyde or ketone groups found besides the carboxylic acid"