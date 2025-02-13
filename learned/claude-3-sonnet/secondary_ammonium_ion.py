"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
neighbors = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in charged_n.GetNeighbors()]