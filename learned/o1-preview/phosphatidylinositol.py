"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
for neighbor in neighbors:
    if neighbor.HasSubstructMatch(inositol_mol):
        phosphate_bonded_to_inositol = True
        break