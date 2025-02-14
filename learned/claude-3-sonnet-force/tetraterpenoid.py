"""
Classifies: CHEBI:26935 tetraterpenoid
"""
largest_chain = max(len(chain) for chain in Chem.Molecules.EnumerateAliphaticSpiroCycles(mol))