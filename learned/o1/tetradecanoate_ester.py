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
    It contains an ester group with a 14-carbon acyl chain derived from tetradecanoic acid.

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
    ester_pattern = Chem.MolFromSmarts('[#6;X3](=O)[O;X2][#6]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"
    
    # For each ester group, check if the acyl chain is tetradecanoyl
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon atom index
        ester_o_idx = match[2]     # Ester oxygen atom index

        carbonyl_c_atom = mol.GetAtomWithIdx(carbonyl_c_idx)
        ester_o_atom = mol.GetAtomWithIdx(ester_o_idx)
        
        # Initialize visited atoms set to prevent infinite loops
        visited_atoms = set()
        visited_atoms.add(ester_o_idx)  # Exclude ester oxygen to avoid traversing into alcohol part
        
        # Function to traverse the acyl chain recursively
        def traverse_acyl_chain(atom, visited_atoms):
            count = 0
            atom_idx = atom.GetIdx()
            visited_atoms.add(atom_idx)
            if atom.GetAtomicNum() == 6:
                count += 1  # Count carbon atoms
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited_atoms and neighbor_idx != ester_o_idx:
                    count += traverse_acyl_chain(neighbor, visited_atoms)
            return count

        # Count carbons in the acyl chain (including carbonyl carbon)
        acyl_chain_length = traverse_acyl_chain(carbonyl_c_atom, visited_atoms)
        
        if acyl_chain_length == 14:
            return True, "Contains tetradecanoate ester group"
    
    return False, "No tetradecanoate ester groups found"