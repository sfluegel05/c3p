"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid has at least 3 carboxyl (-COOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for a terminal carboxyl group
    # - The carbon must have a double bond to O, and a single bond to -OH
    # - The -OH oxygen must be terminal.
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX1H1]")
    if carboxyl_pattern is None:
        return None, "Invalid SMARTS pattern for carboxyl group"

    # Find all matches of the carboxyl group
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    # Check if the number of matches is at least 3
    if len(carboxyl_matches) >= 3:
        return True, "Contains at least three carboxyl groups"
    else:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, less than 3"