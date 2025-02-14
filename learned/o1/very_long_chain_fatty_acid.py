"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: very long-chain fatty acid
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is a fatty acid which has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid group (C(=O)OH or deprotonated form)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[O;H1,-1]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Exclude molecules with ester or amide functional groups
    ester_pattern = Chem.MolFromSmarts('[CX3](=O)[O][CX3,CX4]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if ester_matches:
        return False, "Ester functional group(s) found"

    amide_pattern = Chem.MolFromSmarts('[CX3](=O)[NX3][CX3,CX4]')
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if amide_matches:
        return False, "Amide functional group(s) found"

    # Get the carbon atom of the carboxyl group
    carboxyl_carbons = [match[0] for match in carboxylic_acid_matches]

    # Function to find the longest carbon chain starting from a given atom
    def get_longest_chain(atom_idx, visited):
        max_length = 0
        for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
            neighbor = bond.GetOtherAtomIdx(atom_idx)
            if neighbor not in visited and mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6:
                visited.add(neighbor)
                length = 1 + get_longest_chain(neighbor, visited)
                visited.remove(neighbor)
                if length > max_length:
                    max_length = length
        return max_length

    max_chain_length = 0
    for c_idx in carboxyl_carbons:
        visited = set([c_idx])
        chain_length = get_longest_chain(c_idx, visited)
        if chain_length > max_chain_length:
            max_chain_length = chain_length

    if max_chain_length > 22:
        return True, f"Longest carbon chain length is {max_chain_length}, which is greater than 22"
    else:
        return False, f"Longest carbon chain length is {max_chain_length}, which is not greater than 22"