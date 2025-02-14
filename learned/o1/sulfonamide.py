"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is an amide derivative of a sulfonic acid, characterized by the functional group RS(=O)2NR'2.

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

    # Define sulfonamide SMARTS pattern: sulfur atom double-bonded to two oxygens and single-bonded to nitrogen
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)N")
    if sulfonamide_pattern is None:
        return False, "Invalid SMARTS pattern for sulfonamide"

    # Check if the molecule contains the sulfonamide functional group
    if mol.HasSubstructMatch(sulfonamide_pattern):
        return True, "Contains sulfonamide functional group"
    else:
        return False, "Does not contain sulfonamide functional group"