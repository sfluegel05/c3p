"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: Monoamine, defined as an aralkylamino compound which contains one amino group 
connected to an aromatic ring by a two‐carbon chain.
Examples include tyramine, (R)-noradrenaline, dopamine, etc.
"""

from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as an aralkylamino compound having one amino group connected
    to an aromatic ring by a two-carbon chain. This means there must be a substructure 
    where a non-aromatic nitrogen (which can be neutral or positively charged) is bonded 
    to two saturated carbon atoms (sp3, using [CX4]) that then connect to an aromatic carbon.
    Furthermore, the nitrogen should not be part of an amide bond.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a monoamine, False otherwise.
        str: A message detailing the reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that matches a monoamine-like connectivity:
    #   (1) The nitrogen can be either uncharged or positively charged: 
    #       [$([NX3;!a]),$([N+;!a])]
    #   (2) Followed by two non-aromatic saturated carbons: [CX4][CX4]
    #   (3) Followed by an aromatic carbon: [c]
    # This ensures the detected fragment is: N–C–C–(aromatic C)
    smarts = "[$([NX3;!a]),$([N+;!a])][CX4][CX4][c]"
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        return False, "Error in SMARTS pattern"

    # Find all substructure matches that satisfy the monoamine connectivity.
    matches = mol.GetSubstructMatches(pattern)
    if len(matches) == 0:
        return False, "No aralkylamino moiety (amine connected by a two‐carbon chain to an aromatic ring) found"
    
    # We now filter the matches to exclude those where the nitrogen is part of an amide (or similar)
    # by checking if the nitrogen (first atom in each match) has any neighbor (other than the attached chain)
    # that is in a carbonyl bond (C=O).
    valid_matches = []
    for match in matches:
        n_idx, c1_idx, c2_idx, ar_idx = match
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Assume valid flag; then check all other heavy neighbors of the nitrogen.
        valid = True
        for nb in n_atom.GetNeighbors():
            # Skip the carbon that is part of the detected chain.
            if nb.GetIdx() == c1_idx:
                continue
            # If this neighbor is a carbon that is double-bonded to an oxygen, likely part of an amide,
            # then skip this match.
            if nb.GetAtomicNum() == 6:
                for bond in nb.GetBonds():
                    # Check if bond order is 2 and the neighbor atom is oxygen.
                    # (Use GetDoubleBondEquivalent for robustness.)
                    if bond.GetBeginAtom().GetIdx() == nb.GetIdx():
                        other = bond.GetEndAtom()
                    else:
                        other = bond.GetBeginAtom()
                    if other.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2:
                        valid = False
                        break
            if not valid:
                break
        if valid:
            valid_matches.append(match)
    
    # Now we expect exactly one valid monoamine group.
    if len(valid_matches) == 0:
        return False, "No valid monoamine group detected after filtering out amide-like moieties"
    elif len(valid_matches) > 1:
        return False, f"Found {len(valid_matches)} aralkylamino moieties; expected exactly one monoamine group"
    
    return True, "Contains exactly one amino group connected via a two-carbon chain to an aromatic ring"

# Example test cases (uncomment below to run tests)
# test_smiles = [
#     "CNC[C@@H](O)c1ccc(O)c(O)c1",  # (S)-adrenaline
#     "NCCc1ccc(O)cc1",             # tyramine
#     "[NH3+]CCc1ccc(O)cc1",         # tyraminium (should match monoamine)
#     "OC(=O)[C@H](Cc1ccc(O)c(O)c1)\\N=C/C=C1C[C@H](NC(=C\\1)C(O)=O)C(O)=O"  # dopaxanthin (expected to be rejected)
# ]
# for s in test_smiles:
#     result, reason = is_monoamine(s)
#     print(f"SMILES: {s}\n Result: {result}\n Reason: {reason}\n")