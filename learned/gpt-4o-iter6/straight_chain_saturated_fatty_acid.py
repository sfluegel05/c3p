"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

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
    
    # Check for terminal carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Count carbons and ensure all are sp3 hybridized (except for one carbon in carboxyl group)
    num_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            num_carbons += 1
            # Skip sp3 check for the carboxylic acid carbon
            if (atom.GetDegree() == 3 and any(n.GetAtomicNum() == 8 for n in atom.GetNeighbors())):
                continue
            if atom.GetHybridization() != rdchem.HybridizationType.SP3:
                return False, f"Contains unsaturated carbon (non-sp3 hybridized); SMILES: {smiles}"

    # Ensure the main carbon chain is straight
    longest_chain = Chem.rdmolops.GetSSSR(mol)
    if len(longest_chain) != num_carbons - 1:
        return False, "Chain is not straight"

    return True, f"Contains {num_carbons} carbons in a straight-chain saturated fatty acid structure"