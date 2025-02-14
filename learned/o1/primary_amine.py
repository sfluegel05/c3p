"""
Classifies: CHEBI:32877 primary amine
"""
if neighbor_atom.GetBondBetweenAtoms(neighbor_atom.GetIdx(), atom.GetIdx()).GetBondType() == Chem.BondType.TRIPLE: