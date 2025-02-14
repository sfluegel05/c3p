"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
neighbors = [mol.GetAtomWithIdx(nei).GetAtomicNum() for nei in atom.GetNeighbors()]