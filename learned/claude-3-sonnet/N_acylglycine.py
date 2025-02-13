"""
Classifies: CHEBI:16180 N-acylglycine
"""
common_n_atoms = list(set([atom.GetBeginAtomIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]) & set([atom for match in acyl_atoms for atom in match]) & set([atom for match in glycine_atoms for atom in match]))