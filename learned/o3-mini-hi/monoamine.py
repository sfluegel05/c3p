"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: Monoamine
Definition: A monoamine is an aralkylamino compound having one amino group connected to an aromatic ring
by a two‐carbon chain. In this implementation we search for an aliphatic (non‐aromatic, non‐cyclic) nitrogen,
not involved in an amide bond, that is connected via two consecutive non‐aromatic, acyclic carbons to
an aromatic carbon. The molecule must contain exactly one such aralkylamino moiety.
"""

from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    
    A monoamine is defined as an aralkylamino compound having exactly one amino group connected 
    to an aromatic ring via a two‐carbon chain. In our implementation the algorithm looks for an 
    aliphatic (non-aromatic, non-cyclic) nitrogen not involved in an amide bond that is connected 
    (through two acyclic, non-aromatic carbons) to an aromatic carbon.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a monoamine, False otherwise.
        str: A message explaining the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    monoamine_count = 0  # count how many monoamine moieties are found

    # Iterate over all atoms looking for candidate nitrogen atoms.
    for atom in mol.GetAtoms():
        # Check atom is nitrogen, not aromatic, and not in a ring (i.e. a free amino group)
        if atom.GetAtomicNum() != 7:
            continue
        if atom.GetIsAromatic() or atom.IsInRing():
            continue
        
        # Exclude likely amide-type nitrogens: if any neighboring carbon is doubly bonded to an oxygen,
        # then skip this nitrogen.
        is_amide = False
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 6:  # carbon neighbor
                for bond in nb.GetBonds():
                    # Check if the bond is a double bond.
                    if bond.GetBondTypeAsDouble() == 2:
                        other = bond.GetOtherAtom(nb)
                        if other.GetAtomicNum() == 8:
                            is_amide = True
                            break
                if is_amide:
                    break
        if is_amide:
            continue
        
        # Look for the defined chain: N - C1 - C2 - Ar.
        # We require that both C1 and C2 are aliphatic (i.e. non‐aromatic, non‐cyclic) carbons.
        found_chain = False
        for c1 in atom.GetNeighbors():
            if c1.GetAtomicNum() != 6:  # must be carbon
                continue
            if c1.GetIsAromatic() or c1.IsInRing():
                continue
            # c1 is our first chain carbon.
            for c2 in c1.GetNeighbors():
                if c2.GetIdx() == atom.GetIdx():
                    continue  # avoid backtracking to nitrogen
                if c2.GetAtomicNum() != 6:
                    continue
                if c2.GetIsAromatic() or c2.IsInRing():
                    continue
                # Now look for an aromatic carbon attached to c2.
                for ar in c2.GetNeighbors():
                    if ar.GetIdx() == c1.GetIdx():
                        continue
                    if ar.GetAtomicNum() == 6 and ar.GetIsAromatic():
                        # We found a chain N - C1 - C2 - Ar (exactly three bonds from the candidate N to this aromatic carbon).
                        found_chain = True
                        break
                if found_chain:
                    break
            if found_chain:
                break
                
        if found_chain:
            monoamine_count += 1

    # Exactly one monoamine group should be present.
    if monoamine_count == 0:
        return False, "No aralkylamino moiety (amine connected by a two‐carbon chain to an aromatic ring) found"
    elif monoamine_count > 1:
        return False, f"Found {monoamine_count} aralkylamino moieties; expected exactly one monoamine group"
    else:
        return True, "Contains exactly one amino group connected via a two‐carbon chain to an aromatic ring"
        
# Example usage (for manual testing):
# test_smiles = [
#     "CNC[C@@H](O)c1ccc(O)c(O)c1",  # (S)-adrenaline: valid
#     "NCCc1ccc(O)cc1",              # tyramine: valid
#     "COc1cc(ccc1OCCCN1CCC(CC1)c1noc2cc(F)ccc12)C(C)=O"  # iloperidone: should be accepted if it contains exactly one moiety
# ]
# for s in test_smiles:
#     result, reason = is_monoamine(s)
#     print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")