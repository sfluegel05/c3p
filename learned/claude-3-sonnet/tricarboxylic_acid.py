"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: tricarboxylic acid
Definition: An oxoacid containing three carboxy groups
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid must contain exactly three carboxylic acid (-COOH) groups.

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

    # Look for carboxylic acid groups (-C(=O)OH)
    # Using a SMARTS pattern that specifically looks for -C(=O)OH
    carboxyl_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid groups found"
    
    # Count carboxylic acid groups
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    num_carboxyl = len(carboxyl_matches)
    
    if num_carboxyl != 3:
        return False, f"Found {num_carboxyl} carboxylic acid groups, need exactly 3"

    # Additional validation to ensure we're not counting other groups
    # Count total oxygens and carbons to make sure numbers make sense
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:  # Each COOH has 2 oxygens, so minimum 6 oxygens needed
        return False, "Insufficient oxygen atoms for three carboxylic acids"

    # Make sure the carboxylic groups are not part of esters or anhydrides
    ester_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2][#6]')  # Ester pattern
    anhydride_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2][CX3](=[OX1])')  # Anhydride pattern
    
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    anhydride_matches = len(mol.GetSubstructMatches(anhydride_pattern))
    
    if ester_matches > 0:
        return False, "Contains ester groups"
    if anhydride_matches > 0:
        return False, "Contains acid anhydride groups"

    # Ensure molecule is acidic (has H+ donors)
    if not any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()):
        return False, "No hydrogen atoms found"

    return True, "Contains exactly three carboxylic acid groups"