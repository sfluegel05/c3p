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
    Quinic acid is a cyclitol carboxylic acid with a specific substitution pattern on its cyclohexane core.

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

    # Basic cyclohexane ring with one carboxyl group
    core_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1")
    if core_pattern is None:
        return False, "Invalid SMARTS pattern for core"
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No cyclohexane core found"

    # Check for carboxylic acid or ester group
    carboxyl_pattern = Chem.MolFromSmarts("[$([C](=[O])[OH]),$([C](=[O])[O][C])]")
    if carboxyl_pattern is None:
        return False, "Invalid SMARTS pattern for carboxyl"
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid or ester group found"

    # Count basic features
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if carbon_count < 7:  # minimum carbons in quinic acid
        return False, "Too few carbons for quinic acid structure"
    if oxygen_count < 5:  # minimum oxygens in quinic acid
        return False, "Too few oxygens for quinic acid structure"

    # Look for hydroxyl or ester groups
    hydroxyl_pattern = Chem.MolFromSmarts("[$([OH]),$([O][C](=[O])[#6])]")
    if hydroxyl_pattern is None:
        return False, "Invalid SMARTS pattern for hydroxyl"
    
    # Get the core atoms
    core_matches = mol.GetSubstructMatches(core_pattern)
    if not core_matches:
        return False, "No valid cyclohexane core found"
    
    core_atoms = set(core_matches[0])
    
    # Count oxygen substituents connected to core
    oxygen_substituents = 0
    for atom_idx in core_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                oxygen_substituents += 1

    if oxygen_substituents < 4:  # Need at least 4 oxygen substituents
        return False, "Insufficient oxygen substituents on core"

    # Check for characteristic quinic acid substitution pattern
    # At least one carboxyl and multiple hydroxyls/esters on the ring
    quinic_pattern = Chem.MolFromSmarts("[C]1([C](=[O])[O,H])[C]([O,H])[C]([O,H])[C]([O,H])[C]([O,H])[C]1([O,H])")
    if quinic_pattern is None:
        return False, "Invalid SMARTS pattern for quinic acid"
    if not mol.HasSubstructMatch(quinic_pattern):
        return False, "Does not match quinic acid substitution pattern"

    # Additional check for connected groups
    for atom in mol.GetAtoms():
        # Exclude molecules with phosphate, sulfate, or other non-typical groups
        if atom.GetAtomicNum() not in [1, 6, 8]:  # Only H, C, O allowed
            return False, "Contains non-typical atoms for quinic acid"

    return True, "Contains quinic acid core with characteristic substitution pattern"