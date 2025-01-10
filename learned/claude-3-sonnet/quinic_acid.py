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

    # Look for basic cyclohexane core
    cyclohexane = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1")
    if not mol.HasSubstructMatch(cyclohexane):
        return False, "No cyclohexane ring found"

    # Find cyclohexane rings
    ring_matches = mol.GetSubstructMatches(cyclohexane)
    
    found_valid_core = False
    for ring_match in ring_matches:
        ring_atoms = set(ring_match)
        
        # Check for carboxylic acid or ester attached to ring
        carboxyl_pattern = Chem.MolFromSmarts("[$([CX3](=[OX1])[OX2H1]),$([CX3](=[OX1])[OX2][C])]")
        if not mol.HasSubstructMatch(carboxyl_pattern):
            continue
            
        # Count oxygens attached to ring carbons
        oxygen_count = 0
        for ring_atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(ring_atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    oxygen_count += 1
        
        # Quinic acid should have at least 4 oxygens attached to ring
        # (including the carboxyl group)
        if oxygen_count >= 4:
            found_valid_core = True
            break
    
    if not found_valid_core:
        return False, "No cyclohexane ring with sufficient oxygen substitution found"

    # Additional checks for overall composition
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if carbon_count < 7:  # minimum carbons in quinic acid
        return False, "Too few carbons for quinic acid structure"
    
    if oxygen_count < 5:  # minimum oxygens in quinic acid
        return False, "Too few oxygens for quinic acid structure"

    # Look for ester groups (for derivatives)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_count = len(mol.GetSubstructMatches(ester_pattern))
    
    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    # Must have either free hydroxyls or ester groups
    if hydroxyl_count + ester_count < 3:
        return False, "Insufficient hydroxyl/ester groups"

    return True, "Contains cyclohexane core with carboxylic acid and multiple hydroxyl/ester groups characteristic of quinic acid"