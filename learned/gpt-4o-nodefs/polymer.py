"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    Polymers are characterized by repeating units and long chain structures.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool or None: True if molecule is likely a polymer, False if unlikely, None if classification is not possible.
        str or None: Reason for classification, or None if not classifiable.
    """
    # Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Attempt to find simple chain patterns typical in polymers.
    chain_length_threshold = 25  # Increased threshold for "long chains"
    
    longest_chain = 0
    
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:
            continue  # Skip branch points for simple linear chain searches
        current_chain = 1
        visited_atoms = set()
        
        # Simple DFS to find longest chain
        to_visit = [(atom, 0)]
        while to_visit:
            current, depth = to_visit.pop()
            if depth > current_chain:
                current_chain = depth
            for neighbor in current.GetNeighbors():
                if neighbor.GetIdx() not in visited_atoms:
                    to_visit.append((neighbor, depth + 1))
                    visited_atoms.add(neighbor.GetIdx())
        
        longest_chain = max(longest_chain, current_chain)

    if longest_chain >= chain_length_threshold:
        return True, f"Found long chain (length: {longest_chain}) typical for polymers"

    # If the SMILES string is complex but lacks definitive polymer indicators, return inconclusive.
    atom_counts = mol.GetNumAtoms()
    if atom_counts > 100:
        return None, "Complex structure, inconclusive without further analysis"

    return False, "Does not fit polymer characteristics in basic checks"