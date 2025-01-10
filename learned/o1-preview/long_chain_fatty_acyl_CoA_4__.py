"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
coa_pattern = Chem.MolFromSmarts("NC(=O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](CO[P](=O)(O)OP(=O)(O)OP(=O)(O)OC[C@H]2O[C@H](n3cnc4c(N)ncnc43)[C@H](O)[C@H]2O)[C@H](O)[C@@H]1O")