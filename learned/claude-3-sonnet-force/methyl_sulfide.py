"""
Classifies: CHEBI:86315 methyl sulfide
"""
neighbors = [mol.GetAtomWithIdx(neighbor_idx) for neighbor_idx in sulfur.GetNeighbors()]