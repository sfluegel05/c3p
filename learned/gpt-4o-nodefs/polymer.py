"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    Polymers are characterized by repeating units (monomers) and specific structural features.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is likely a polymer, False if unlikely.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define potential polymer SMARTS patterns
    repeating_unit_patterns = [
        Chem.MolFromSmarts("C(C(C)=O)"),  # ester unit -COOH (e.g., acrylates)
        Chem.MolFromSmarts("CNC(=O)"),   # amide unit -CONH- (e.g., nylons)
        Chem.MolFromSmarts("COC"),       # ether unit -O- (e.g., polyethylene oxide)
    ]

    # Search for repeating units
    for pattern in repeating_unit_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches polymer-related structure: {Chem.MolToSmarts(pattern)}"

    # Check chain length
    chain_length_threshold = 30  # Dynamic threshold based on polymer features
    longest_chain = 0
    
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:
            continue  # Skip branch points for simple linear chain searches
        current_chain = 1
        visited_atoms = set()
        
        # Simple DFS to find longest chain
        to_visit = [(atom, 1)]
        while to_visit:
            current, depth = to_visit.pop()
            visited_atoms.add(current.GetIdx())
            if depth > current_chain:
                current_chain = depth
            for neighbor in current.GetNeighbors():
                if neighbor.GetIdx() not in visited_atoms:
                    to_visit.append((neighbor, depth + 1))
                    visited_atoms.add(neighbor.GetIdx())
        
        longest_chain = max(longest_chain, current_chain)

    if longest_chain >= chain_length_threshold:
        return True, f"Found long chain (length: {longest_chain}) typical for polymers"

    return False, "Does not fit polymer characteristics"