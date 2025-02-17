"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
#!/usr/bin/env python
"""
Classifies: Nucleoside phosphate
A nucleobase‐containing molecular entity that is a nucleoside in which 
one or more of the sugar hydroxy groups has been converted into a mono‐ or poly‐phosphate.
This includes both nucleotides and “non‐nucleotide” nucleoside phosphates.

Strategy revision:
  1. Parse the SMILES and require at least one phosphorus atom.
  2. Find a candidate sugar: a non‐aromatic 5‐membered ring containing 1 oxygen and 4 carbons.
  3. Check for a nucleobase attached to the sugar. Instead of only checking immediate neighbors,
     we also allow a one-atom (e.g. CH₂) linker. We require that the “base” is part of an aromatic ring system
     (even if fused) containing at least 2 nitrogen atoms.
  4. Verify that a phosphate is attached to the sugar in one of the following ways:
       – a sugar atom is directly bonded to a phosphorus, or
       – a non–ring oxygen substituent of a sugar atom is directly bound to a phosphorus, or
       – a bridging (exocyclic) carbon (e.g. the 5′–CH₂ group) carries an oxygen attached to phosphorus.
       
If all conditions are met, we return True.
Note: There are many edge cases so this is one possible implementation.
"""
from rdkit import Chem

def is_nucleobase_candidate(mol, atom_idx, visited=None, depth=0, max_depth=1):
    """
    Recursively check from a starting atom (typically a neighbor of the sugar)
    if we can reach an atom that is part of an aromatic ring (of typical size 5- or 6-membered, 
    or fused system) which contains at least 2 nitrogen atoms. 
    We allow traversing one extra bond if needed (max_depth=1).
    """
    if visited is None:
        visited = set()
    # Prevent revisiting atoms
    if atom_idx in visited:
        return False
    visited.add(atom_idx)
    
    atom = mol.GetAtomWithIdx(atom_idx)
    ring_info = mol.GetRingInfo()
    # Check every ring that contains this atom:
    for ring in ring_info.AtomRings():
        # We relax the strict size check to allow fused rings.
        if len(ring) < 5 or len(ring) > 10:
            continue
        # Check if all atoms in this ring are aromatic.
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetIsAromatic() for a in ring_atoms):
            continue
        # Count nitrogen atoms in the ring.
        n_nitrogen = sum(1 for a in ring_atoms if a.GetAtomicNum() == 7)
        if n_nitrogen >= 2:
            return True
    # If not found and we have not exceeded depth, try one more bond.
    if depth < max_depth:
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in visited:
                if is_nucleobase_candidate(mol, nbr.GetIdx(), visited, depth+1, max_depth):
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
    
    # 1. Check that the molecule has at least one phosphorus atom.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_atoms:
        return False, "No phosphorus atoms found; unlikely to be a nucleoside phosphate"
    
    ring_info = mol.GetRingInfo()
    sugar_ring = None
    # 2. Find candidate sugar ring: non–aromatic 5–membered ring with 1 oxygen and 4 carbons.
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(atom.GetIsAromatic() for atom in atoms):
            continue
        num_oxygens = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        num_carbons = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if num_oxygens == 1 and num_carbons == 4:
            sugar_ring = ring
            break
    if sugar_ring is None:
        return False, "No suitable sugar ring (non‐aromatic 5‐membered ring with 1 O and 4 C) found"
    
    sugar_set = set(sugar_ring)
    
    # 3. Check for a nucleobase attached to the sugar.
    nucleobase_found = False
    # In typical nucleosides the anomeric carbon carries the nucleobase.
    # Here, scan every atom in the sugar ring.
    for idx in sugar_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # If neighbor is in the sugar, skip.
            if nbr.GetIdx() in sugar_set:
                continue
            # Check if the neighbor (or one bond further) is part of an aromatic ring with ≥2 N.
            if is_nucleobase_candidate(mol, nbr.GetIdx()):
                nucleobase_found = True
                break
        if nucleobase_found:
            break
    if not nucleobase_found:
        return False, "No nucleobase found attached to the sugar (no aromatic heterocycle with ≥2 N found nearby)"
    
    # 4. Check for a phosphate group attached to the sugar.
    phosphate_attached = False
    for idx in sugar_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # Skip atoms that belong to the sugar ring.
            if nbr.GetIdx() in sugar_set:
                continue
            # Case A: Direct bond to phosphorus.
            if nbr.GetAtomicNum() == 15:
                phosphate_attached = True
                break
            # Case B: An oxygen substituent directly bound to phosphorus.
            if nbr.GetAtomicNum() == 8:
                for oxy_nbr in nbr.GetNeighbors():
                    if oxy_nbr.GetAtomicNum() == 15:
                        phosphate_attached = True
                        break
                if phosphate_attached:
                    break
            # Case C: A bridging carbon (exocyclic, e.g. CH2) with an oxygen that is bound to phosphorus.
            if nbr.GetAtomicNum() == 6:
                for sub in nbr.GetNeighbors():
                    # Skip if it is the original sugar atom or part of the sugar ring.
                    if sub.GetIdx() == idx or sub.GetIdx() in sugar_set:
                        continue
                    if sub.GetAtomicNum() == 8:
                        for oxy_nbr in sub.GetNeighbors():
                            if oxy_nbr.GetAtomicNum() == 15:
                                phosphate_attached = True
                                break
                        if phosphate_attached:
                            break
                if phosphate_attached:
                    break
        if phosphate_attached:
            break
    if not phosphate_attached:
        return False, "No phosphate group attached to the sugar (expected via oxygen or bridging carbon)"
    
    return True, "Contains a furanose sugar with an attached nucleobase and a phosphate (direct or via a CH2 group) indicative of a nucleoside phosphate"

# Example usage:
if __name__ == "__main__":
    # Test one of the given molecules: 2'-deoxy-5-methyl-5'-cytidylic acid
    test_smiles = "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N"
    result, reason = is_nucleoside_phosphate(test_smiles)
    print(result, reason)