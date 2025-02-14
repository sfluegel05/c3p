"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: CHEBI:82198 wax
"""

from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is defined as an organic compound or mixture of compounds
    that is composed of long-chain molecules and is malleable at ambient temperatures.
    Commonly, waxes are esters formed from long-chain fatty acids and long-chain alcohols.

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

    # Define ester pattern
    ester_pattern = Chem.MolFromSmarts("[C](=O)O[C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester groups found"

    # Define threshold for long-chain (number of carbons)
    chain_length_threshold = 12

    # Check chain lengths for each ester group
    for match in ester_matches:
        # Atom indices for the ester group
        C_carbonyl_idx = match[0]  # Carbonyl carbon
        O_ester_idx = match[2]     # Ester oxygen
        C_alcohol_idx = match[3]   # Carbon attached to ester oxygen

        # Get chain length starting from carbonyl carbon
        acid_chain_length = get_chain_length(mol, C_carbonyl_idx, O_ester_idx)
        # Get chain length starting from alcohol carbon
        alcohol_chain_length = get_chain_length(mol, C_alcohol_idx, O_ester_idx)

        if acid_chain_length >= chain_length_threshold and alcohol_chain_length >= chain_length_threshold:
            return True, f"Ester with long chains found (acid chain: {acid_chain_length} carbons, alcohol chain: {alcohol_chain_length} carbons)"

    return False, "No ester with long chains found"

def get_chain_length(mol, start_atom_idx, exclude_atom_idx):
    """
    Counts the number of carbons in the chain starting from start_atom_idx, 
    excluding the bond to exclude_atom_idx. Traverses until chain ends or branches.

    Args:
        mol: RDKit molecule object
        start_atom_idx: Starting atom index
        exclude_atom_idx: Atom index to exclude (previous atom)

    Returns:
        int: Number of carbons in the chain
    """
    chain_length = 0
    current_atom_idx = start_atom_idx
    prev_atom_idx = exclude_atom_idx
    while True:
        atom = mol.GetAtomWithIdx(current_atom_idx)
        if atom.GetAtomicNum() != 6:
            break  # Stop if not carbon
        chain_length += 1
        # Get neighbors excluding previous atom
        neighbors = [n for n in atom.GetNeighbors() if n.GetIdx() != prev_atom_idx]
        if len(neighbors) != 1:
            # Termination (either end of chain or branching)
            break
        prev_atom_idx = current_atom_idx
        current_atom_idx = neighbors[0].GetIdx()
    return chain_length