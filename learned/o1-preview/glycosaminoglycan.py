"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
sugar_smarts = Chem.MolFromSmarts("""
    [C;R]1([O;R][C;R][C;R][C;R][C;R]1) | 
    [C;R]1([O;R][C;R][C;R][C;R]1)
""")