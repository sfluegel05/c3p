"""
Classifies: CHEBI:87657 octanoate ester
"""
# Exclude ester oxygen and any heteroatoms (to ensure unbranched carbon chain)
for nb in neighbors:
    nb_idx = nb.GetIdx()
    nb_atom_num = nb.GetAtomicNum()

    if nb_idx == ester_o_idx:
        continue  # Skip ester oxygen
    elif nb_atom_num == 6:
        stack.append((nb_idx, current_atom_idx))
    else:
        return False, "Acyl chain is substituted or contains heteroatoms"