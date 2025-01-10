"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: CHEBI:36586 2-oxo monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid is any monocarboxylic acid having a 2-oxo substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"

    # Check for a 2-oxo substituent (carbonyl group within two bonds of the carboxylic acid)
    # The pattern looks for a carbonyl group (C=O) within two bonds of the carboxylic acid
    oxo_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4][CX4][CX3](=O)[OX2H1]")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if len(oxo_matches) == 0:
        # If no matches, try a more flexible pattern that allows for a single bond between the carbonyl and carboxylic acid
        oxo_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4][CX3](=O)[OX2H1]")
        oxo_matches = mol.GetSubstructMatches(oxo_pattern)
        if len(oxo_matches) == 0:
            return False, "No 2-oxo substituent found"

    return True, "Contains a single carboxylic acid group and a 2-oxo substituent"