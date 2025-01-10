"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: CHEBI:16865 quinic acid and its derivatives
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or a quinic acid derivative based on its SMILES string.
    Quinic acid is a cyclitol carboxylic acid with a cyclohexane core containing multiple hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid derivative, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclohexane core with carboxylic acid
    cyclohexane_carboxylic = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1[C](=O)[OH]")
    if not mol.HasSubstructMatch(cyclohexane_carboxylic):
        return False, "No cyclohexane ring with carboxylic acid group found"

    # Count hydroxyl groups (both free and esterified)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    ester_count = len(mol.GetSubstructMatches(ester_pattern))
    
    total_oxy_groups = hydroxyl_count + ester_count
    
    # Quinic acid derivatives should have at least 3 oxygen-containing groups
    # (including both free hydroxyls and ester groups) attached to the cyclohexane ring
    if total_oxy_groups < 3:
        return False, f"Insufficient oxygen-containing groups (found {total_oxy_groups}, need at least 3)"

    # Check for basic carbon framework
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 7:  # minimum carbons in quinic acid
        return False, "Too few carbons for quinic acid structure"

    # Additional check for cyclohexane ring with multiple oxygens
    cyclohexane_with_oxygens = Chem.MolFromSmarts("[C]1[C]([O])[C]([O])[C]([O])[C]([O])[C]1")
    if not mol.HasSubstructMatch(cyclohexane_with_oxygens):
        return False, "Cyclohexane ring lacks required oxygen substitution pattern"

    # Check if it's a simple carboxylic acid (must have COOH group)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern) and not mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=[OX1])[OX2][C]")):
        return False, "No carboxylic acid or ester group found"

    return True, "Contains cyclohexane core with carboxylic acid and multiple hydroxyl/ester groups characteristic of quinic acid"