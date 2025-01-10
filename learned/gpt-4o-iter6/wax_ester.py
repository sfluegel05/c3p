"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of the 
    carboxy group of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.

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

    # Ester linkage pattern: C(=O)O
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    # Find ester group and analyze the carbon chains extending from it
    matches = mol.GetSubstructMatches(ester_pattern)

    for match in matches:
        carbonyl_carbon = match[0]
        ester_oxygen = match[2]
        
        side1_chain_length = _get_chain_length(mol, carbonyl_carbon)
        side2_chain_length = _get_chain_length(mol, ester_oxygen)
        
        # Consider a total minimum chain length threshold
        total_chain_length = side1_chain_length + side2_chain_length

        # Relax the minimum requirement for individual chains, while emphasizing total
        if total_chain_length >= 20 and side1_chain_length >= 8 and side2_chain_length >= 8:
            return True, "Molecule contains a fatty acid ester linkage with adequate carbon chain length"
    
    return False, "Carbon chains are too short to be considered fatty acid/alcohol"

def _get_chain_length(mol, start_atom_index):
    """
    Recursively finds the length of a carbon chain starting from a given atom.
    Only considers carbon atoms in a linear sequence.
    """
    visited = set()
    max_chain_length = 0

    def _traverse_chain(atom_index):
        nonlocal max_chain_length
        if atom_index in visited:
            return 0
        visited.add(atom_index)

        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetAtomicNum() != 6:  # Only consider carbon atoms
            return 0

        local_chain_length = 1  # Count the current carbon atom
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited:
                local_chain_length += _traverse_chain(neighbor.GetIdx())
        
        max_chain_length = max(max_chain_length, local_chain_length)
        return local_chain_length

    _traverse_chain(start_atom_index)
    return max_chain_length