"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    Ultra-long-chain fatty acids have a chain length greater than C27 and a carboxylic acid group.

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

    # Look for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Find the longest carbon chain
    num_carbons = 0
    for chain in Chem.rdmolops.GetSymmSSSR(mol):
        atom_indices = list(chain)
        carbon_count = sum(1 for i in atom_indices if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
        num_carbons = max(num_carbons, carbon_count)

    # Check if chain length is greater than C27
    if num_carbons > 27:
        return True, f"Contains carboxylic acid group and chain length is C{num_carbons}, which is greater than C27"
    else:
        return False, f"Chain length is C{num_carbons}, which is not greater than C27"