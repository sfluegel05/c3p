"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is a sugar alcohol with multiple hydroxyl groups and no carbonyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of carbonyl groups (C=O) or similar structures
    carbonyl_pattern = Chem.MolFromSmarts("[$([CX3]=[OX1]),$([CX3+]=[OX1-])]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains carbonyl-like group(s), not an alditol"

    # Count hydroxyl groups (O connected to any C and also to H) which should be plentiful in alditols
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need at least 3 to be an alditol"

    # Ensure the molecule primarily consists of carbon and hydroxyl groups
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8, 1):  # Carbon, Oxygen, Hydrogen
            return False, f"Contains atoms other than C, O, and H: {atom.GetSymbol()}"

    # Check that molecule doesn't possess other unexpected functions
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings > 0:
        return False, "Alditol should not contain rings; potentially cyclic sugars detected"

    return True, "Contains multiple hydroxyl groups and lacks carbonyl groups, consistent with alditol structure."