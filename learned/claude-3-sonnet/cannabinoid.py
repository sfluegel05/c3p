"""
Classifies: CHEBI:67194 cannabinoid
"""
aromatic_rings = [ring for ring in AllChem.GetSymmSSSR(mol) if ring.IsAromatic()]