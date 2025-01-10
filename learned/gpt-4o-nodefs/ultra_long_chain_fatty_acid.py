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

    # Count the carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 28:
        return False, f"Not enough carbon atoms for ultra-long-chain fatty acid: {c_count}"

    # Look for long carbon chain
    chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")  # At least 9 CH2 in sequence
    num_long_chains = len(mol.GetSubstructMatches(chain_pattern))
    
    if num_long_chains < 2:  # Preferably two or more segments of long chains, considering branching or other motifs
        return False, f"Insufficient long hydrocarbon chains: {num_long_chains}"

    # Check for characteristic ULCFA substructures
    polyunsaturated_pattern = Chem.MolFromSmarts("C=C")
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    methoxy_pattern = Chem.MolFromSmarts("[C][O][C]")
    cyclopropyl_pattern = Chem.MolFromSmarts("C1CC1")

    has_required_groups = any(mol.HasSubstructMatch(pattern) for pattern in [polyunsaturated_pattern, hydroxyl_pattern, methoxy_pattern, cyclopropyl_pattern])

    if has_required_groups:
        return True, "Contains ultra-long-chain fatty acid features"

    return False, "Insufficient structural features for an ultra-long-chain fatty acid"