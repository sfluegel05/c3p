"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
path = Chem.rdmolops.GetShortestPath(mol, list(isoquinoline_atoms[i])[0], list(isoquinoline_atoms[j])[0])