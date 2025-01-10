"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
if rdMolDescriptors.CalcMolFormula(mol).count('C') < 28:
    return False, "Insufficient carbon content for ultra-long chain fatty acid"