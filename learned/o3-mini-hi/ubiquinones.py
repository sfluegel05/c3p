"""
Classifies: CHEBI:16389 ubiquinones
"""
#!/usr/bin/env python
"""
Classifies: Ubiquinones – natural benzoquinones derived from 2,3-dimethoxy-5-methylbenzoquinone.
They share a redox-active quinone ring (a cyclohexadiene-1,4-dione core) that typically carries at least one methoxy
substituent and a methyl group on the ring. If any long substituent chain (≥5 carbons) is attached to the core,
it must contain at least one C=C double bond (i.e. be a polyprenyl chain) rather than a simple alkyl chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is classified as a ubiquinone.
    A ubiquinone must have:
      - A quinone core: a six-membered ring that contains exactly 2 carbonyl (C=O) groups.
      - The two carbonyl groups should be at “nonadjacent” (ideally para) positions in the ring.
      - At least one methoxy substituent (an oxygen directly attached to the ring that itself is bound to a methyl)
        on one of the ring atoms.
      - At least one methyl (CH3) group directly attached to the core.
      - If any substituent chain (other than small groups) is attached to the core and contains 5 or more carbons,
        the chain must contain a C=C double bond (to qualify as a polyprenyl side chain).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True and an explanation if classified as a ubiquinone, otherwise False and a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information and look for a six-membered ring that might be the quinone core.
    rings = mol.GetRingInfo().AtomRings()
    core_ring = None
    core_idx_set = None
    for ring in rings:
        if len(ring) != 6:
            continue
        # Check how many ring atoms are "carbonyl carbons" connected via a double bond to an oxygen.
        carbonyl_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                # Check for a double bond (C=O) with an oxygen that is not in the ring.
                if bond.GetBondTypeAsDouble() == 1:
                    other = bond.GetOtherAtom(atom)
                    if other.GetSymbol() == "O" and other.GetIdx() not in ring:
                        carbonyl_indices.append(idx)
                        break
        # Accept the ring as core if it has exactly 2 carbonyl groups.
        if len(carbonyl_indices) == 2:
            # As a simple check, ensure the two carbonyl carbons are non-adjacent in the ring.
            # In a 6-membered ring, para positions are separated by 3 atoms.
            # We check the ring order (sorted order is not helpful so we use index order as in the ring list)
            pos = [ring.index(ci) for ci in carbonyl_indices]
            pos.sort()
            # Calculate cyclic distance
            dist = min(pos[1]-pos[0], 6 - (pos[1]-pos[0]))
            if dist >= 3:
                core_ring = ring
                core_idx_set = set(ring)
                break
    if core_ring is None:
        return False, "Quinone core (six-membered ring with two carbonyl groups in para positions) not found"

    # Initialize counters for substituents.
    methoxy_count = 0  # groups like -OCH3 attached directly.
    methyl_count = 0   # CH3 groups directly attached.
    identified_substituents = set()  # atom indices that are already classified as small substituents.

    # Helper function: check if an atom qualifies as a CH3 group.
    def is_methyl(atom):
        # Atom is carbon and has at least 3 implicit or explicit hydrogens.
        return (atom.GetSymbol() == "C" and atom.GetTotalNumHs() >= 3)
    
    # Examine each atom in the core ring.
    for idx in core_idx_set:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in core_idx_set:
                continue
            # Check for methoxy: neighbor is oxygen.
            if nbr.GetSymbol() == "O":
                nbr_neighbors = [a for a in nbr.GetNeighbors() if a.GetIdx() != idx]
                if len(nbr_neighbors) == 1:
                    attached = nbr_neighbors[0]
                    if is_methyl(attached):
                        methoxy_count += 1
                        identified_substituents.add(nbr_idx)
                        identified_substituents.add(attached.GetIdx())
            # If neighbor is carbon and not already identified, check for methyl.
            elif nbr.GetSymbol() == "C":
                if nbr_idx in identified_substituents:
                    continue
                if is_methyl(nbr):
                    methyl_count += 1
                    identified_substituents.add(nbr_idx)

    # Now check for any substituent chain attached to the core which wasn't already classified as a small group.
    polyprenyl_found = False
    polyprenyl_reason = ""
    visited_substitution = set()
    chain_error = None

    def dfs_chain(start_atom, visited):
        # Depth-first-search to get the substituent fragment.
        stack = [start_atom]
        frag_atoms = set()
        has_double = False
        while stack:
            a = stack.pop()
            if a.GetIdx() in visited:
                continue
            visited.add(a.GetIdx())
            frag_atoms.add(a.GetIdx())
            for bond in a.GetBonds():
                nb = bond.GetOtherAtom(a)
                if nb.GetIdx() in core_idx_set:
                    continue
                if nb.GetIdx() not in visited:
                    stack.append(nb)
                if bond.GetBondTypeAsDouble() == 1:
                    has_double = True
        return frag_atoms, has_double

    # Iterate over neighbors of core atoms (exclude those already identified as methoxy or methyl)
    for idx in core_idx_set:
        core_atom = mol.GetAtomWithIdx(idx)
        for nbr in core_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in core_idx_set or nbr_idx in identified_substituents or nbr_idx in visited_substitution:
                continue
            # Only consider carbon chains.
            if nbr.GetSymbol() != "C":
                continue
            frag, has_dbl = dfs_chain(nbr, visited_substitution)
            n_carbons = sum(1 for i in frag if mol.GetAtomWithIdx(i).GetSymbol() == "C")
            if n_carbons >= 5:
                # For a long chain we require at least one double bond.
                if has_dbl:
                    polyprenyl_found = True
                    polyprenyl_reason = f"Found a substituent chain with {n_carbons} carbons containing C=C double bonds (indicative of a polyprenyl chain)."
                else:
                    chain_error = f"Found a long substituent chain with {n_carbons} carbons but with no C=C double bonds (likely a non-isoprenoid alkyl chain)."
                break
        if chain_error or polyprenyl_found:
            break

    # Check the minimal requirements.
    if methoxy_count < 1:
        return False, f"Found {methoxy_count} methoxy group(s) on the quinone core; need at least 1."
    if methyl_count < 1:
        return False, "No methyl group directly attached to the quinone core found."
    if chain_error:
        return False, chain_error

    # Build a summary explanation.
    reason = f"Quinone core detected (six-membered ring with two carbonyls in para positions). Found {methoxy_count} methoxy group(s) and {methyl_count} methyl group(s). "
    if polyprenyl_found:
        reason += polyprenyl_reason
    else:
        reason += "No qualifying long polyprenyl side chain detected (this is acceptable for some ubiquinone homologues)."

    return True, reason

# Example usage for testing:
if __name__ == "__main__":
    test_smiles = [
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-2
        "C1(=C(C(C(=C(C1=O)OC)[O-])=O)C)C/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)\\C)\\C)\\C)\\C",  # 3-demethyl-ubiquinone-7(1-)
        "COC1=C(OC)C(=O)C(C)=CC1=O",  # ubiquinone-0
        "COC1=C(OC)C(=O)C(CC=C(C)C)=C(C)C1=O",  # ubiquinone-1
        "COC1=C(OC)C(=O)C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)=C(C)C1=O",  # ubiquinone-8
    ]
    
    for s in test_smiles:
        valid, explanation = is_ubiquinones(s)
        print(f"SMILES: {s}\nResult: {valid}\nExplanation: {explanation}\n")