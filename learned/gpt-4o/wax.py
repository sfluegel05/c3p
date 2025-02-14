"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is considered a wax based on its SMILES string.
    A wax is characterized by long-chain molecules, typically containing esters,
    being malleable at ambient temperatures, and having significant flexibility.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a wax, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for ester linkage (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    # Initiate count arrays for carbon linkages
    c_counts = [0] * mol.GetNumAtoms()
    max_chain_length = 0
    
    # A utility function to perform depth-first search for longest carbon chain
    def dfs_longest_chain(atom_index, current_length, visited):
        nonlocal max_chain_length
        visited[atom_index] = True
        for neighbor in mol.GetAtomWithIdx(atom_index).GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 6 and not visited[neighbor_idx]:
                dfs_longest_chain(neighbor_idx, current_length + 1, visited)
        max_chain_length = max(max_chain_length, current_length)
        visited[atom_index] = False
    
    # Check all carbon atoms for longest chain path
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            visited = [False] * mol.GetNumAtoms()
            dfs_longest_chain(atom.GetIdx(), 1, visited)
    
    if max_chain_length < 14:
        return False, f"Longest carbon chain is {max_chain_length}, which is too short for a wax"
    
    # Calculate number of rotatable bonds to ascertain flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, f"Insufficient rotatable bonds ({n_rotatable}), waxes should be flexible"
    
    return True, f"Structure meets criteria for wax: Contains ester linkage, long carbon chain ({max_chain_length} carbons), and sufficient rotatable bonds ({n_rotatable})"