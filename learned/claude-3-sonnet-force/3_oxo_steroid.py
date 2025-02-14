"""
Classifies: CHEBI:47788 3-oxo steroid
"""
oxo_neighbor_bond = mol.GetBondBetweenAtoms(oxo_atom_idx, oxo_neighbor_idx)
if oxo_neighbor_bond.GetBondDir() != BondDir.ENDDOWNRIGHT: