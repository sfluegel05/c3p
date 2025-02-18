"""
Classifies: CHEBI:24026 fatty alcohol
"""
def get_longest_chain(mol):
    # Get all carbon atoms
    carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    max_length = 0
    for start in carbons:
        for end in carbons:
            if start == end:
                continue
            path = Chem.GetShortestPath(mol, start, end)
            if len(path) > max_length:
                max_length = len(path)
    return max_length