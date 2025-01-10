"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Count number of carbons and ensure all are sp3 hybridized (indicating single bonds)
    num_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            num_carbons += 1
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                return False, f"Contains unsaturated carbon (non-sp3 hybridized); SMILES: {smiles}"
    
    # Check for straight chain: Each carbon should have 2 hydrogen atoms (except functional groups)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            num_bonds = sum(1 for _ in atom.GetBonds())
            if num_bonds > 2 and not any(n.GetAtomicNum() == 8 for n in atom.GetNeighbors()):  # Carbon not in carboxylic group
                return False, "Carbon atom forms more than 2 bonds with non-terminal atoms"

    return True, f"Contains {num_carbons} carbons in a straight-chain saturated fatty acid structure"