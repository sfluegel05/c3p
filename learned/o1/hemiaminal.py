"""
Classifies: CHEBI:73080 hemiaminal
"""
hemiaminal_pattern = Chem.MolFromSmarts('[C;X4;!$(C=O)](-[OH])(-[N;!$(N-C=[O,N]);!$(N=[C,N,O]);!a])')