"""
Classifies: CHEBI:73754 thiosugar
"""
neighbors = [mol.GetAtomWithIdx(nbr_idx) for nbr_idx in sulfur.GetNeighbors()]