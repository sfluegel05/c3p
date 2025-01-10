"""
Classifies: CHEBI:32878 alkene
"""
total_bonds = sum(bond.GetBondTypeAsDouble() for bond in mol.GetBonds())
expected_bonds = (2 * c_count - 1)  # For alkene: single bonds + one double bond