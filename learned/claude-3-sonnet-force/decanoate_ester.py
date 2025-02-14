"""
Classifies: CHEBI:87658 decanoate ester
"""
if any(mol.GetBondBetweenAtoms(decanoate_idx[1], alcohol_idx[0]).GetBondType() == Chem.BondType.ESTER):