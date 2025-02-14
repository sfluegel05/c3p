"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: CHEBI:24062 sulfonamide
An amide of a sulfonic acid RS(=O)2NR'2, where R' can be H, alkyl, aryl, etc.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.

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
    
    # Look for sulfonamide groups
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)N[!H]")  # Sulfonamide pattern, excluding NH groups
    
    sulfonamide_matches = mol.GetSubstructMatches(sulfonamide_pattern)
    
    if not sulfonamide_matches:
        return False, "No sulfonamide groups found"
    
    return True, "Contains one or more sulfonamide groups"