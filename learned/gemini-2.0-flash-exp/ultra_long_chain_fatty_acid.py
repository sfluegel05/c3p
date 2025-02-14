"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid
    (chain length > C27) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Molecule does not contain a carboxylic acid group."

    # Find the carboxyl carbon
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    carboxyl_carbon_idx = matches[0][0] # Get the index of the carbon in C=O
    
    def get_longest_chain(start_atom_idx, current_path, max_path_len):
        max_path_len[0] = max(max_path_len[0], len(current_path))
        
        atom = mol.GetAtomWithIdx(start_atom_idx)
        
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 6 and neighbor_idx not in current_path:
                get_longest_chain(neighbor_idx, current_path + [neighbor_idx], max_path_len)
        return
    
    max_path_len = [0]
    get_longest_chain(carboxyl_carbon_idx, [carboxyl_carbon_idx], max_path_len)
    
    # Count carbons in longest path
    longest_chain_carbon_count = max_path_len[0]

    if longest_chain_carbon_count > 27:
        return True, f"Longest carbon chain has {longest_chain_carbon_count} carbons, which is > 27."
    else:
        return False, f"Longest carbon chain has {longest_chain_carbon_count} carbons, which is <= 27."