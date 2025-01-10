"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    Polymers are characterized by repeating units (monomers), specific structural features,
    and often have long chain-like structures with possible branching or network formation.

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

    # Define a list of potential repeating unit SMARTS patterns
    repeating_unit_patterns = [
        Chem.MolFromSmarts("C(C(C)=O)"),  # ester unit -COOH (e.g., acrylates)
        Chem.MolFromSmarts("CNC(=O)"),   # amide unit -CONH- (e.g., nylons),
        Chem.MolFromSmarts("COC"),       # ether unit -O- (e.g., polyethylene oxide)
        Chem.MolFromSmarts("C(=C)C"),    # vinyl unit (e.g., polyvinyl chloride)
        Chem.MolFromSmarts("[*]([*])[*]") # generic branching pattern
    ]

    # Check for repeating units indicative of polymer structures
    for pattern in repeating_unit_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches polymer-related structure: {Chem.MolToSmarts(pattern)}"

    # Calculate longest chain and check if it's indicative of a polymer
    chain_length_threshold = 30  # Adjusted threshold
    longest_chain = 0
    
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:
            continue  # Skip obvious branch points for simple linear chains
        visited_atoms = set()
        to_visit = [(atom, 0)]  # Include branching depth

        while to_visit:
            current, depth = to_visit.pop()
            if current.GetIdx() in visited_atoms:
                continue
            visited_atoms.add(current.GetIdx())
            longest_chain = max(longest_chain, depth)
            for neighbor in current.GetNeighbors():
                if neighbor.GetIdx() not in visited_atoms:
                    if neighbor.GetDegree() > 2:
                        branch_depth = depth + 1
                    else:
                        branch_depth = depth + 1
                    to_visit.append((neighbor, branch_depth))

    if longest_chain >= chain_length_threshold:
        return True, f"Found long chain (length: {longest_chain}) typical for polymers"

    return False, "Does not fit polymer characteristics"