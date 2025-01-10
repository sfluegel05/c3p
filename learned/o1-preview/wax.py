"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    Waxes are esters formed from long-chain fatty acids and long-chain alcohols,
    typically with carbon chains of 14 carbons or more on both sides of the ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify ester functional groups with atom mapping
    ester_pattern = Chem.MolFromSmarts("[C:1](=O)[O:2][C:3]")  # Ester pattern with atom mapping
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional group found"

    for match in ester_matches:
        # Get atom indices
        carbonyl_c_idx = match[0]
        single_o_idx = match[1]
        alkyl_c_idx = match[2]

        # Function to count carbons on acyl side (from carbonyl carbon)
        def count_acyl_carbons(atom_idx, exclude_idx=None, visited=None):
            if visited is None:
                visited = set()
            if atom_idx in visited:
                return 0
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            count = 0
            if atom.GetAtomicNum() == 6:  # Carbon atom
                count += 1
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Exclude the oxygen atom connected in the ester group
                if neighbor_idx != exclude_idx and neighbor_idx != single_o_idx:
                    count += count_acyl_carbons(neighbor_idx, atom_idx, visited)
            return count

        # Function to count carbons on alkoxy side (from alkyl carbon attached to oxygen)
        def count_alkoxy_carbons(atom_idx, exclude_idx=None, visited=None):
            if visited is None:
                visited = set()
            if atom_idx in visited:
                return 0
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            count = 0
            if atom.GetAtomicNum() == 6:  # Carbon atom
                count += 1
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Exclude the oxygen atom connected in the ester group
                if neighbor_idx != exclude_idx and neighbor_idx != single_o_idx:
                    count += count_alkoxy_carbons(neighbor_idx, atom_idx, visited)
            return count

        # Count carbons on acyl side
        acyl_chain_length = count_acyl_carbons(carbonyl_c_idx, single_o_idx)
        # Count carbons on alkoxy side
        alkoxy_chain_length = count_alkoxy_carbons(alkyl_c_idx, single_o_idx)

        if acyl_chain_length >= 14 and alkoxy_chain_length >= 14:
            return True, f"Ester with acyl chain length {acyl_chain_length} and alkoxy chain length {alkoxy_chain_length}"

    return False, "Ester present, but chains are not long enough for wax"

__metadata__ = {
    'chemical_class': {
        'name': 'wax',
        'definition': 'A chemical substance that is an organic compound or mixture of compounds that is composed of long-chain molecules and is malleable at ambient temperatures.'
    }
}