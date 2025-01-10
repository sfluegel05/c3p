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
    if c_count < 30:
        return False, f"Not enough carbon atoms for ultra-long-chain fatty acid: {c_count}"

    # Look for long carbon chain and verify bonds
    chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")  # Simplified long chain matcher
    has_long_chain = mol.HasSubstructMatch(chain_pattern)
    
    if not has_long_chain:
        return False, "Lacks a sufficiently long hydrocarbon chain"

    # Look for additional functional groups or structural motifs
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    methoxy_pattern = Chem.MolFromSmarts("[C][O][C]")
    cyclopropyl_pattern = Chem.MolFromSmarts("C1CC1")
    
    has_important_groups = any(mol.HasSubstructMatch(pattern) for pattern in [hydroxyl_pattern, methoxy_pattern, cyclopropyl_pattern])
    
    if has_important_groups:
        return True, "Contains ultra-long-chain fatty acid features"

    return False, "Insufficient structural features for an ultra-long-chain fatty acid"