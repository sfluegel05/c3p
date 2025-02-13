"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is considered a wax based on its SMILES string.
    A wax is characterized by long-chain hydrocarbons, typically containing esters,
    being malleable at ambient temperatures, favoring linear aliphatic chains.
    
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
    
    # Finding and evaluating the longest carbon chains linked to ester linkage
    carbon_chains = []
    
    # Find ester group atoms
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    for est_match in ester_matches:
        explored = set()
        def explore_chain(atom_idx, length=0):
            explored.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                length += 1
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx not in explored and neighbor.GetAtomicNum() == 6:
                    length = explore_chain(n_idx, length)
            return length
        
        # Check both sides of the ester linkage (C and O atom indices)
        chain_lengths = [explore_chain(est_match[0]), explore_chain(est_match[2])]
        carbon_chains.extend(chain_lengths)
    
    # Check for at least one long chain
    if any(chain >= 16 for chain in carbon_chains):
        return True, "Contains ester linkage and at least one long carbon chain indicative of a wax"
    
    return False, "No sufficiently long carbon chains found for wax"