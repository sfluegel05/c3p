"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: wax
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    Waxes are esters formed from long-chain fatty acids and long-chain alcohols, typically with
    carbon chains of 14 carbons or more.

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

    # Find ester functional groups
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional group found"

    # For each ester group, check the length of the carbon chains on both sides
    for match in ester_matches:
        # Get the carbonyl carbon atom index
        c_atom_idx = match[0]
        # Get the oxygen atom index of the ester bond
        o_atom_idx = match[2]
        
        # Get the acyl chain (fatty acid side)
        acyl_chain_length = get_chain_length(mol, c_atom_idx, exclude_atom_idx=o_atom_idx)
        # Get the alkoxy chain (fatty alcohol side)
        alkoxy_chain_length = get_chain_length(mol, o_atom_idx, exclude_atom_idx=c_atom_idx)
        
        # Check if both chains are long (e.g., at least 14 carbons)
        if acyl_chain_length >= 14 and alkoxy_chain_length >= 14:
            return True, f"Ester with acyl chain length {acyl_chain_length} and alkoxy chain length {alkoxy_chain_length}"
    
    return False, "Ester present, but chains are not long enough for wax"

def get_chain_length(mol, start_atom_idx, exclude_atom_idx):
    """
    Calculates the length of the carbon chain starting from the given atom index,
    excluding the direction towards the exclude_atom_idx.

    Args:
        mol (Chem.Mol): RDKit molecule object
        start_atom_idx (int): Atom index to start the chain from
        exclude_atom_idx (int): Atom index to exclude from traversal (e.g., the ester bond)

    Returns:
        int: Number of carbons in the chain
    """
    visited = set()
    chain_length = 0
    stack = [(start_atom_idx, None)]
    while stack:
        current_atom_idx, previous_atom_idx = stack.pop()
        if current_atom_idx in visited or current_atom_idx == exclude_atom_idx:
            continue
        visited.add(current_atom_idx)
        atom = mol.GetAtomWithIdx(current_atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_length += 1
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx != previous_atom_idx:
                    stack.append((neighbor_idx, current_atom_idx))
    return chain_length

__metadata__ = {
    'chemical_class': {
        'name': 'wax',
        'definition': 'A chemical substance that is an organic compound or mixture of compounds that is composed of long-chain molecules and is malleable at ambient temperatures.'
    }
}