"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
backbone = mol.GetSubstructMatch(Chem.MolFragmentToSmarts(mol, atomsToUse=path))