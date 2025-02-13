"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:26356 diol
A compound that contains two hydroxy groups, generally assumed to be, but not necessarily, alcoholic.
Aliphatic diols are also called glycols.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    n_hydroxy = len(hydroxy_matches)
    
    # Diols must have exactly 2 hydroxy groups
    if n_hydroxy != 2:
        return False, f"Found {n_hydroxy} hydroxy groups, need exactly 2"
    
    # Check for hydroxy groups in specific positions
    # e.g., cyclic structures, vicinal/geminal diols, etc.
    # Using SMARTS patterns or other substructure matching techniques
    
    # Example: Vicinal diol pattern
    vicinal_diol_pattern = Chem.MolFromSmarts("[OX2H][CX4H2][OX2H]")
    if mol.HasSubstructMatch(vicinal_diol_pattern):
        return True, "Contains vicinal diol substructure"
    
    # Add more specific patterns as needed
    
    # Check for restricted functional groups (if desired)
    # e.g., aldehydes, ketones, etc.
    
    # Example: Exclude aldehydes
    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
    if mol.HasSubstructMatch(aldehyde_pattern):
        return False, "Molecule contains aldehyde group"
    
    # If no specific patterns match, return a general classification
    return True, "Contains two hydroxy groups"