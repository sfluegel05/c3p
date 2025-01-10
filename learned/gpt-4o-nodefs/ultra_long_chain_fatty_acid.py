"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a carboxylic acid group (C(O)=O)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 22:
        return False, f"Not enough carbon atoms for ultra-long-chain fatty acid: {c_count}"

    # Look for continuous carbon chain; match ~C6 or longer (covers both straight and with minimal branching or unsaturation)
    chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    long_chain_matches = mol.GetSubstructMatches(chain_pattern)
    
    total_long_chain_atoms = sum(len(match) for match in long_chain_matches)
    if total_long_chain_atoms < 16:  # Ensure overall contribution to a lengthy span
        return False, f"Insufficient long hydrocarbon chain contribution: {total_long_chain_atoms}"
    
    # Check for characteristic ULCFA substructures
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")  # OH group
    methoxy_pattern = Chem.MolFromSmarts("[C][O][C]")  # Methoxy group
    cyclopropyl_pattern = Chem.MolFromSmarts("C1CC1")  # Cyclopropyl ring

    has_functional_groups = any(mol.HasSubstructMatch(pattern) for pattern in [
        hydroxyl_pattern, methoxy_pattern, cyclopropyl_pattern
    ])
    
    if has_functional_groups:
        return True, "Contains ultra-long-chain fatty acid features"

    return False, "Insufficient structural features for an ultra-long-chain fatty acid"