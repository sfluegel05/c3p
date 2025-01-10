"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is characterized by a terminal aldehyde attached to a long
    aliphatic carbon chain, which may have double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to get the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to identify a terminal aldehyde group: alkyl chain ending with =O
    # It looks for a carbonyl carbon attached to another carbon
    aldehyde_pattern = Chem.MolFromSmarts('O=[C;!R]') # Ensures non-branching, non-ring
    
    # Check the main features
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if not aldehyde_matches:
        return False, "No terminal aldehyde group found"
    
    # Check for qualifiable aliphatic chains, looking for long chain structures
    for match in aldehyde_matches:
        aldehyde_atom_idx = match[1] # The carbon atom in [C]=O

        # Traverse from aldehyde C to assess chain structure
        def is_valid_aliphatic(atom):
            # Carbon atoms that are not in rings and are part of the main chain
            return atom.GetAtomicNum() == 6 and not atom.IsInRing()

        visited = set()
        chain_length = 0
        
        def traverse_chain(atom_idx, prev_idx):
            nonlocal chain_length
            if atom_idx not in visited:
                visited.add(atom_idx)
                atom = mol.GetAtomWithIdx(atom_idx)
                if is_valid_aliphatic(atom):
                    chain_length += 1
                    # Traverse neighbors except the one we came from
                    for neighbor in atom.GetNeighbors():
                        next_idx = neighbor.GetIdx()
                        if next_idx != prev_idx:
                            traverse_chain(next_idx, atom_idx)

        # Start from the connected carbon, and ensure a long chain
        traverse_chain(aldehyde_atom_idx - 1, -1)  # Go from the non-oxygen side

        if chain_length >= 6:  # Typical chain length criteria for fatty structures
            return True, "Valid fatty aldehyde: Terminal aldehyde group with a suitable aliphatic chain"
    
    return False, "Carbon chain too short for typical fatty aldehyde"