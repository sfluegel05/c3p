"""
Classifies: CHEBI:35366 fatty acid
"""
longest_chain_length = max(len(set([mol.GetAtomWithIdx(idx).GetNeighbors()[0].GetIdx()
                                    for match in chain_matches
                                    for idx in match])) + 1, default=0)