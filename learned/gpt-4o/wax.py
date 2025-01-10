"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is typically an ester formed from long-chain fatty acids and long-chain alcohols, with one or two ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester groups (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester group found"
    if len(ester_matches) > 2:
        return False, f"Too many ester groups, found {len(ester_matches)}"

    # Identify all carbon atoms excluding ester groups
    visited_atoms = set()
    all_chain_lengths = []

    def compute_chain_length(atom, visited_atoms):
        # Recursive function that computes chain length, avoiding atoms in ester groups
        visited_atoms.add(atom.GetIdx())
        chain_lengths = []
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited_atoms:
                if not mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[O-1]")):  # Check for ester linkage
                    chain_length = compute_chain_length(neighbor, visited_atoms) + 1
                    chain_lengths.append(chain_length)
        return max(chain_lengths, default=0)

    # Calculate maximum chain length for each carbon
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() not in visited_atoms:
            chain_length = compute_chain_length(atom, set())
            all_chain_lengths.append(chain_length)

    # Two longest carbon chains need to be significant in length (e.g., more than 14 carbons)
    long_chain_threshold = 14
    significant_long_chains = sorted(all_chain_lengths, reverse=True)[:2]
    
    if len(significant_long_chains) < 2:
        return False, "Insufficient long carbon chains; need at least two"
    if any(chain < long_chain_threshold for chain in significant_long_chains):
        return False, "Longest chains are too short; need at least two chains over 14 carbons"

    return True, "Contains long-chain molecules with one or two ester linkages characteristic of waxes"