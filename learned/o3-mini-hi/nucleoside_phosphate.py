"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
#!/usr/bin/env python
"""
Classifies: Nucleoside phosphate
A nucleobase‐containing molecular entity that is a nucleoside in which 
one or more of the sugar hydroxy groups has been converted into a mono‐ or poly‐phosphate.
This includes both nucleotides and “non‐nucleotide” nucleoside phosphates.

Revised strategy:
  1. Parse the SMILES and require at least one phosphorus atom.
  2. Find a candidate sugar ring: we require a non‐aromatic 5‐membered ring with exactly one oxygen and four carbons.
  3. Identify an anomeric carbon in the sugar ring that has an external substituent.
     That substituent (even via one CH₂ linker) should lead to an aromatic heterocycle containing at least 2 nitrogen atoms.
  4. Verify that a phosphate group (P atom) is attached to the sugar ring either directly,
     or via a non‐ring oxygen substituent or via a bridging carbon bearing an –OH that in turn bonds to P.
If all conditions are met, we return True.
Note: There remain many edge cases so this is one possible implementation.
"""
from rdkit import Chem

def is_nucleobase_candidate(mol, start_idx, max_depth=1):
    """
    Starting from start_idx (which is attached to the sugar) see if within max_depth
    we can reach a fully aromatic ring that contains at least 2 nitrogen atoms.
    This function uses a simple breadth-first search.
    """
    seen = {start_idx}
    current_level = [start_idx]
    for depth in range(max_depth+1):
        next_level = []
        for idx in current_level:
            atom = mol.GetAtomWithIdx(idx)
            # Look at each ring that contains this atom.
            for ring in mol.GetRingInfo().AtomRings():
                # We allow rings with 5-10 atoms (to cover fused systems)
                if len(ring) < 5 or len(ring) > 10:
                    continue
                # Only consider rings that are entirely aromatic.
                if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    continue
                # Count the number of nitrogen atoms.
                n_nitrogen = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)
                if n_nitrogen >= 2:
                    return True
            # Expand search
            for nbr in atom.GetNeighbors():
                nxt = nbr.GetIdx()
                if nxt not in seen:
                    seen.add(nxt)
                    next_level.append(nxt)
        if next_level:
            current_level = next_level
        else:
            break
    return False


def phosphate_linked_to_sugar(mol, sugar_ring):
    """
    Check that at least one atom of the sugar ring has a substituent that is either:
      A. A phosphorus atom;
      B. An oxygen atom directly connected to a phosphorus;
      C. A bridging carbon (likely CH2) with an oxygen that is connected to a phosphorus.
    """
    sugar_set = set(sugar_ring)
    for idx in sugar_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in sugar_set:
                continue
            # Condition A: Direct phosphorus connection
            if nbr.GetAtomicNum() == 15:
                return True
            # Condition B: Oxygen substituent that in turn connects to phosphorus.
            if nbr.GetAtomicNum() == 8:
                for oxy_nbr in nbr.GetNeighbors():
                    if oxy_nbr.GetIdx() == idx or oxy_nbr.GetIdx() in sugar_set:
                        continue
                    if oxy_nbr.GetAtomicNum() == 15:
                        return True
            # Condition C: Bridging carbon possibility.
            if nbr.GetAtomicNum() == 6:
                for sub in nbr.GetNeighbors():
                    if sub.GetIdx() in (idx, nbr.GetIdx()) or sub.GetIdx() in sugar_set:
                        continue
                    if sub.GetAtomicNum() == 8:
                        for sub_nbr in sub.GetNeighbors():
                            if sub_nbr.GetIdx() in (nbr.GetIdx(), sub.GetIdx()):
                                continue
                            if sub_nbr.GetAtomicNum() == 15:
                                return True
    return False


def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    Args:
       smiles (str): SMILES string of the molecule.
    Returns:
       (bool, str): Tuple (True, reason) if the molecule meets classification criteria,
                    else (False, explanation).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # 1. Require at least one phosphorus atom.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_atoms:
        return False, "No phosphorus atoms found; unlikely to be a nucleoside phosphate"
        
    ring_info = mol.GetRingInfo()
    candidate_sugar = None
    # 2. Find candidate sugar ring: non‐aromatic, 5‐membered ring with exactly one O and four C.
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(a.GetIsAromatic() for a in atoms):
            continue
        num_oxygens = sum(1 for a in atoms if a.GetAtomicNum() == 8)
        num_carbons = sum(1 for a in atoms if a.GetAtomicNum() == 6)
        if num_oxygens == 1 and num_carbons == 4:
            candidate_sugar = ring
            break
    if candidate_sugar is None:
        return False, "No suitable sugar ring (non‐aromatic 5‐membered ring with 1 O and 4 C) found"
    
    # 3. Identify a connection from the sugar ring to a nucleobase.
    sugar_set = set(candidate_sugar)
    nucleobase_found = False
    # Look for an "anomeric" sugar carbon: an sp3 carbon in the ring with at least one neighbor
    # outside the sugar that (directly or via one extra atom) belongs to an aromatic ring with ≥2 N.
    for idx in candidate_sugar:
        atom = mol.GetAtomWithIdx(idx)
        # Only consider carbon atoms in the sugar (this is typical for the anomeric carbon).
        if atom.GetAtomicNum() != 6:
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in sugar_set:
                continue
            # Check if the neighbor itself is the nucleobase candidate.
            if is_nucleobase_candidate(mol, nbr.GetIdx(), max_depth=0):
                nucleobase_found = True
                break
            # Alternatively, if the neighbor is a CH2 (or similar, i.e. carbon not aromatic),
            # check one more bond away.
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                for sub in nbr.GetNeighbors():
                    if sub.GetIdx() in sugar_set or sub.GetIdx() == atom.GetIdx():
                        continue
                    if is_nucleobase_candidate(mol, sub.GetIdx(), max_depth=0):
                        nucleobase_found = True
                        break
            if nucleobase_found:
                break
        if nucleobase_found:
            break
    if not nucleobase_found:
        return False, "No nucleobase found attached to the sugar (no aromatic heterocycle with ≥2 N found nearby)"
        
    # 4. Check that a phosphate group is attached to the sugar.
    if not phosphate_linked_to_sugar(mol, candidate_sugar):
        return False, "No phosphate group attached to the sugar (expected via oxygen or bridging carbon)"
    
    return True, "Contains a furanose sugar with an attached nucleobase and a phosphate (direct or via a CH2 group) indicative of a nucleoside phosphate"


# Example usage:
if __name__ == "__main__":
    # Test with one of the provided molecules: 2'-deoxy-5-methyl-5'-cytidylic acid
    test_smiles = "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N"
    result, reason = is_nucleoside_phosphate(test_smiles)
    print(result, reason)