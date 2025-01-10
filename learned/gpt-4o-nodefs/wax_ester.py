"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    Wax esters consist of long-chain fatty acids and fatty alcohols 
    connected through ester linkages, typically featuring a significant ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ester linkage pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # If no ester groups are found, it cannot be a wax ester
    if not ester_matches:
        return False, "No ester linkage found"

    # Analyze each ester linkage
    for match in ester_matches:
        c_atom, o1_atom, o2_atom = match

        # Look at the carbon chains attached to both sides of the ester linkage
        carbon_chain_length_1 = 1  # Including the ester carbon
        carbon_chain_length_2 = 1  # Including the ester oxygen for the second chain

        visited = {c_atom, o1_atom, o2_atom}

        # Explore carbon chain 1 (from ester carbon)
        for bond in mol.GetAtomWithIdx(c_atom).GetBonds():
            if bond.GetOtherAtomIdx(c_atom) not in visited:
                carbon_chain_length_1 += expand_chain_length(mol, bond.GetOtherAtomIdx(c_atom), visited)

        # Explore carbon chain 2 (from ester oxygen, traverse the atom attached to oxygen)
        for bond in mol.GetAtomWithIdx(o2_atom).GetBonds():
            if bond.GetOtherAtomIdx(o2_atom) not in visited:
                carbon_chain_length_2 += expand_chain_length(mol, bond.GetOtherAtomIdx(o2_atom), visited)

        # A potential wax ester needs substantial chains on both sides of the ester group
        if carbon_chain_length_1 >= 10 and carbon_chain_length_2 >= 10:
            return True, "Contains ester linkage with sufficient long carbon chains"

    return False, "Ester linkages found, but chains are not sufficiently long"

def expand_chain_length(mol, atom_idx, visited):
    """Recursive function to calculate the chain length from a starting atom."""
    length = 0
    if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6:  # Only count carbons
        length += 1

    visited.add(atom_idx)

    for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
        next_atom = bond.GetOtherAtomIdx(atom_idx)
        if next_atom not in visited:
            length += expand_chain_length(mol, next_atom, visited)

    return length