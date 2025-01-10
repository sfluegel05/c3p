"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
if not Chem.rdmolops.GetMolFrags(mol, sanitized=False, asMols=False, frags=...)[1] == 1:
    return False, "Molecule is not a single continuous structure"