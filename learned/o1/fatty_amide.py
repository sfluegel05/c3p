"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI:24038 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define amide group pattern where nitrogen is not in a ring
    amide_pattern = Chem.MolFromSmarts("C(=O)N([H])[^R]")
    
    # Find all amide groups in the molecule
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No primary amide functional group found"

    # Loop through amide groups to find fatty acyl chains
    for amide_match in amide_matches:
        carbonyl_c_idx = amide_match[0]
        n_idx = amide_match[2]

        # Exclude amide nitrogens attached to rings
        n_atom = mol.GetAtomWithIdx(n_idx)
        if n_atom.IsInRing():
            continue

        # Exclude molecules where nitrogen has more than one substituent (secondary or tertiary amides)
        if n_atom.GetDegree() > 1:
            continue

        # Function to recursively traverse the acyl chain
        def traverse_chain(atom_idx, visited):
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom_idx in visited:
                return 0
            if atom.IsInRing():
                return 0  # Exclude rings
            visited.add(atom_idx)
            if atom.GetAtomicNum() != 6:
                return 0  # Only follow carbons
            length = 1  # Count current carbon
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() not in visited and nbr.GetIdx() != carbonyl_c_idx]
            for nbr_idx in neighbors:
                length += traverse_chain(nbr_idx, visited)
                break  # Only follow linear path (avoid branching)
            return length

        # Find the carbon attached to the carbonyl carbon (start of acyl chain)
        carbonyl_c_atom = mol.GetAtomWithIdx(carbonyl_c_idx)
        acyl_chain_length = 0
        for neighbor in carbonyl_c_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx != n_idx and neighbor.GetAtomicNum() == 6:
                visited_atoms = set()
                acyl_chain_length = traverse_chain(neighbor_idx, visited_atoms)
                break  # Only consider the chain directly attached to the carbonyl carbon

        if acyl_chain_length >= 4:
            return True, f"Molecule is a fatty amide with acyl chain length {acyl_chain_length}"
        else:
            return False, f"Acyl chain too short ({acyl_chain_length} carbons); not a fatty amide"

    return False, "No fatty acyl amide group found"

# Example usage:
# smiles = "CCCCCCCC(=O)NCCO"  # N-decanoylglycine
# result, reason = is_fatty_amide(smiles)
# print(result, reason)