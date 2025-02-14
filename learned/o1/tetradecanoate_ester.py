"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: CHEBI:35980 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is formed by the condensation of tetradecanoic acid with an alcohol or phenol.
    It contains a linear, unbranched acyl chain of 14 carbons (including the carbonyl carbon) connected via an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define ester functional group pattern
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"
    
    # For each ester group, check if the acyl chain is tetradecanoyl
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon atom index
        ester_o_idx = match[1]     # Ester oxygen atom index
        o_adjacent_c_idx = match[2]  # Carbon attached to ester oxygen (alcohol part)

        carbonyl_c_atom = mol.GetAtomWithIdx(carbonyl_c_idx)
        ester_o_atom = mol.GetAtomWithIdx(ester_o_idx)
        
        # Initialize visited atoms set to prevent infinite loops
        visited_atoms = set()
        visited_atoms.add(ester_o_idx)  # Exclude ester oxygen to avoid traversing into alcohol part
        visited_atoms.add(o_adjacent_c_idx)  # Exclude alcohol part

        # Function to traverse the acyl chain linearly
        def traverse_acyl_chain(atom, prev_atom_idx, chain_length):
            atom_idx = atom.GetIdx()
            visited_atoms.add(atom_idx)
            
            # Check if atom is a carbon
            if atom.GetAtomicNum() != 6:
                return None  # Non-carbon atom encountered
            
            # Increment chain length
            chain_length += 1

            # Get neighbor carbons excluding previous atom and ester oxygen
            neighbor_carbons = []
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx != prev_atom_idx and neighbor.GetAtomicNum() == 6:
                    neighbor_carbons.append(neighbor)
                elif neighbor_idx == ester_o_idx:
                    continue  # Ignore ester oxygen
                elif neighbor_idx in visited_atoms:
                    return None  # Cycle detected
            
            # Check for branching
            if len(neighbor_carbons) > 1:
                return None  # Branching occurred
            
            elif len(neighbor_carbons) == 0:
                # Reached terminal carbon
                return chain_length

            else:
                # Continue traversal
                next_atom = neighbor_carbons[0]
                return traverse_acyl_chain(next_atom, atom_idx, chain_length)
        
        # Start traversal from carbonyl carbon, with previous atom index as ester oxygen
        acyl_chain_length = traverse_acyl_chain(carbonyl_c_atom, ester_o_idx, 0)
        
        if acyl_chain_length == 14:
            return True, "Contains tetradecanoate ester group"
    
    return False, "No tetradecanoate ester groups found"