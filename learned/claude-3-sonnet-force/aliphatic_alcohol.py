"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
sum(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in atom.GetNeighbors())