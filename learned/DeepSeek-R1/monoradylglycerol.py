"""
Classifies: CHEBI:76575 monoradylglycerol
"""
glycerol_patterns = [
    Chem.MolFromSmarts("[C](-[OH])(-[OH])-[C]-[C]-O-*"),  # substituent on third carbon
    Chem.MolFromSmarts("[C](-[OH])-[C](-O-*)-[C](-[OH])"),  # substituent on second carbon
    Chem.MolFromSmarts("[C](-O-*)-[C](-[OH])-[C](-[OH])")  # substituent on first carbon
]