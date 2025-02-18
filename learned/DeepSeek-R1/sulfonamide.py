"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is an amide of sulfonic acid with structure RS(=O)₂NR'₂.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define sulfonamide pattern: S connected to two O (double bonds) and one N
    # Accounts for various substitutions on the nitrogen
    sulfonamide_pattern = Chem.MolFromSmarts('[S](=[O])(=[O])-[#7]')
    
    # Check for presence of sulfonamide group
    if mol.HasSubstructMatch(sulfonamide_pattern):
        return True, "Contains sulfonamide group (S(=O)(=O)-N)"
    else:
        return False, "No sulfonamide group detected"