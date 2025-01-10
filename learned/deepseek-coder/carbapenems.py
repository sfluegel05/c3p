"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: CHEBI:60853 carbapenem
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    A carbapenem is a beta-lactam antibiotic with a carbapenem skeleton, which is a bicyclic structure
    containing a beta-lactam ring fused to a five-membered ring, and is substituted at positions 3, 4, and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core carbapenem skeleton pattern
    carbapenem_pattern = Chem.MolFromSmarts("[C@@H]12[C@@H](C)C(=O)N1C=C[C@@H]2C")
    if not mol.HasSubstructMatch(carbapenem_pattern):
        return False, "No carbapenem skeleton found"

    # Check for common substituents at positions 3, 4, and 6
    # Position 3: Typically has a sulfur-containing group (e.g., thioether)
    sulfur_pattern = Chem.MolFromSmarts("[SX2]")
    sulfur_matches = mol.GetSubstructMatches(sulfur_pattern)
    if len(sulfur_matches) == 0:
        return False, "No sulfur-containing substituent at position 3"

    # Position 4: Typically has a hydroxyl group or similar
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl group at position 4"

    # Position 6: Typically has a carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) == 0:
        return False, "No carboxyl group at position 6"

    # Additional checks for beta-lactam ring
    beta_lactam_pattern = Chem.MolFromSmarts("[C@@H]1[C@@H](C)C(=O)N1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    return True, "Contains carbapenem skeleton with typical substituents at positions 3, 4, and 6"