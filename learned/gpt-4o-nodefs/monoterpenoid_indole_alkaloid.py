"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    Monoterpenoid indole alkaloids typically contain an indole moiety and a monoterpenoid structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS patterns for common structural features
    indole_pattern = Chem.MolFromSmarts("c1cc2c([nH]c1)cccc2")  # Indole ring
    monoterpenoid_pattern = Chem.MolFromSmarts("C1CC2CCCC1C2")  # Monoterpenoid unit (cyclized C9 or C10)

    # Check for indole moiety
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole moiety found"

    # Check for monoterpenoid skeleton
    if not mol.HasSubstructMatch(monoterpenoid_pattern):
        return False, "No monoterpenoid skeleton found"

    return True, "Contains indole moiety and monoterpenoid skeleton"