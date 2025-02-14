"""
Classifies: CHEBI:16385 organic sulfide
"""
sulfide_pattern = Chem.MolFromSmarts('[#16D2;$([#16]([#6;!$(C=[O,N,S,P])])[#6;!$(C=[O,N,S,P])]);!$([#16]-[!#6])]')