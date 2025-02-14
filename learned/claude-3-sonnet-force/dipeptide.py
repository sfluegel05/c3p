"""
Classifies: CHEBI:46761 dipeptide
"""
has_carboxyl = any(atom.GetSymbol() == "O" and sum(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in mol.GetAtomNeighbors(atom.GetIdx())) > 0 for atom in mol.GetAtoms())
has_amino = any(atom.GetSymbol() == "N" and sum(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in mol.GetAtomNeighbors(atom.GetIdx())) > 1 for atom in mol.GetAtoms())