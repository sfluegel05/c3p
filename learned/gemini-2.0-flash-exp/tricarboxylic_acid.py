"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid has exactly 3 carboxyl (-COOH) groups.

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

    # Define SMARTS for carboxyl group, including possible protonation or metal forms
    # C(=O)[OH1,OX2-,OM] - where OH1 is the protonated carboxyl group, OX2- is the deprotonated group, and OM is the metal carboxylate.
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH1,OX2-,OM]")

    # Find all matches of the carboxyl group
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    # Check if the number of matches is exactly 3
    if len(carboxyl_matches) == 3:
        return True, "Contains exactly three carboxyl groups"
    else:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, not 3"