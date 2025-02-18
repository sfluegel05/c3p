"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid (CHEBI:32312)
Any straight-chain, C18 polyunsaturated fatty acid having two C=C double bonds.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a straight-chain C18 fatty acid with two double bonds and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for carboxylic acid group
    carb_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carb_acid_pattern):
        return False, "No carboxylic acid group"

    # Total carbons must be 18
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 18:
        return False, f"Expected 18 carbons, found {c_count}"

    # Find the carbonyl carbon of the carboxylic acid
    matches = mol.GetSubstructMatches(carb_acid_pattern)
    if not matches:
        return False, "No carboxylic acid group found"
    carbonyl_atom = mol.GetAtomWithIdx(matches[0][0])

    # Function to find the longest chain starting from the carbonyl atom
    def find_longest_chain(start_atom):
        visited = set()
        longest_path = []

        def dfs(current, path):
            nonlocal longest_path
            if current.GetIdx() in visited:
                return
            visited.add(current.GetIdx())
            path.append(current)
            if len(path) > len(longest_path):
                longest_path = list(path)
            for neighbor in current.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    dfs(neighbor, path)
            path.pop()
            visited.remove(current.GetIdx())

        dfs(start_atom, [])
        return longest_path

    longest_chain = find_longest_chain(carbonyl_atom)
    if len(longest_chain) < 18:
        return False, f"Longest chain from carbonyl has {len(longest_chain)} carbons, need 18"

    # Check chain continuity and count double bonds
    double_bonds_in_chain = 0
    for i in range(len(longest_chain) - 1):
        bond = mol.GetBondBetweenAtoms(longest_chain[i].GetIdx(), longest_chain[i+1].GetIdx())
        if not bond:
            return False, "Discontinuous chain"
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bonds_in_chain += 1

    if double_bonds_in_chain != 2:
        return False, f"Found {double_bonds_in_chain} double bonds in main chain"

    # Check for branching in the main chain
    for i in range(1, len(longest_chain)-1):
        atom = longest_chain[i]
        neighbors_in_chain = 0
        for neighbor in atom.GetNeighbors():
            if neighbor in longest_chain:
                neighbors_in_chain += 1
        if neighbors_in_chain != 2:
            return False, "Branching in main chain"

    return True, "Straight-chain C18 fatty acid with two double bonds and a carboxylic acid group"