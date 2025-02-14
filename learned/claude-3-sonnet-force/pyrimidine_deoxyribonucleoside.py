"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
# Look for pyrimidine rings
pyrimidine_rings = [ring for ring in AllChem.GetSymmSSSR(mol) if sum(1 for atom in ring if atom.GetAtomicNum() in [7, 6]) == 5]
if not pyrimidine_rings:
    return False, "No pyrimidine ring found"