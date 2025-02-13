"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
if atom.GetSymbol() == "C" and atom.GetIsInRingSize() > 6:
    continue  # Ignore large ring carbons