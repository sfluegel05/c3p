"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is characterized by a series of repeating units (monomers).

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

    # Polymers usually have repeating units. However, straightforward identification may not be possible in practice.
    # Check for some indications of repetition. E.g., repeated patterns not easily inferred from SMILES.
    
    # Attempt to find simple patterns or flag potential inconclusive cases.
    # Strategy might be domain-specific based on known patterns; here is a simple placeholder example.

    # Number of rings/cyclic structures can provide a (very rough) hint.
    num_rings = mol.GetRingInfo().NumRings()
    if num_rings > 5:
        return False, "Multiple rings found, likely not a simple polymer"
    
    # Check for long chain structures, which might hint at polymer-like properties
    chain_length_threshold = 10  # Arbitrary threshold for "long chains"
    atom_counts = mol.GetNumAtoms()
    longest_chain = 0
    
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:  # Branch points usually not contributing to linear polymer chains
            continue
        current_chain = 1
        visited_atoms = set()
        
        # Perform simple exploration to find chain length from current atom
        to_visit = [(atom, 0)]
        while to_visit:
            current, depth = to_visit.pop()
            if depth > current_chain:
                current_chain = depth
            for neighbor in current.GetNeighbors():
                if neighbor.GetIdx() not in visited_atoms:
                    to_visit.append((neighbor, depth + 1))
                    visited_atoms.add(neighbor.GetIdx())
        
        if current_chain > chain_length_threshold:
            longest_chain = max(longest_chain, current_chain)

    if longest_chain >= chain_length_threshold:
        return True, f"Found long chain (length: {longest_chain}) typical for polymers"

    # If the SMILES string is highly complex and doesn't fit simple checks, return inconclusive.
    if atom_counts > 50:
        return None, "Complex structure, inconclusive without further analysis"

    return False, "Does not fit polymer characteristics in basic checks"