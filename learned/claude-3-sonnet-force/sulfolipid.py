"""
Classifies: CHEBI:61384 sulfolipid
"""
lipid_chain_pattern = Chem.MolFromSmarts(f"[CX4]~[CX4]~[CX4]~[CX4]~{[carbon_atom.GetIdx()]}")