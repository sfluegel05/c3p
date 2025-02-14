"""
Classifies: CHEBI:50523 butenolide
"""
if ring_info.AtomRings()[0].Count() != 5:
    return False, "Ring is not a 5-membered gammma-lactone"