"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
#!/usr/bin/env python
"""
Classifies: Nucleoside phosphate
A nucleobase-containing molecular entity that is a nucleoside in which 
one or more of the sugar hydroxy groups has been converted into a mono- or poly‐phosphate.
The term includes both nucleotides and non‐nucleotide nucleoside phosphates.
Improved strategy:
  1. Parse the SMILES.
  2. Ensure at least one phosphorus atom is present.
  3. Identify a candidate sugar ring: a non‐aromatic 5‐membered ring with exactly one oxygen and four carbons.
  4. Confirm that a nucleobase is attached to the sugar: a neighboring ring (5- or 6-membered) that is aromatic and has at least 2 nitrogen atoms.
  5. Verify that a phosphate group is attached to the sugar by checking for an exocyclic oxygen (a neighbor of the sugar ring atom not in the ring) that is bonded to a phosphorus.
"""

from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple of classification boolean and an explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the presence of phosphorus atoms.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_atoms:
        return False, "No phosphorus atoms found; unlikely to be a nucleoside phosphate"
    
    # 2. Identify candidate sugar ring.
    # We look for a non-aromatic 5-membered ring with exactly 1 oxygen and 4 carbons.
    ring_info = mol.GetRingInfo()
    sugar_ring = None
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(atom.GetIsAromatic() for atom in atoms):
            continue  # skip if aromatic
        num_oxygens = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        num_carbons = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if num_oxygens == 1 and num_carbons == 4:
            sugar_ring = ring
            break
    if sugar_ring is None:
        return False, "No suitable sugar ring (non-aromatic 5-membered ring with 1 oxygen and 4 carbons) found"
    
    # 3. Identify a nucleobase attached to the sugar.
    # Search neighbors of sugar atoms that are part of an aromatic 5- or 6-membered ring containing ≥2 nitrogen atoms.
    nucleobase_found = False
    sugar_set = set(sugar_ring)
    for idx in sugar_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in sugar_set:
                continue  # skip atoms within the sugar ring
            # Check each ring the neighbor participates in.
            for ring in ring_info.AtomRings():
                if nbr_idx not in ring or len(ring) not in (5, 6):
                    continue
                ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                if not all(a.GetIsAromatic() for a in ring_atoms):
                    continue
                n_nitrogen = sum(1 for a in ring_atoms if a.GetAtomicNum() == 7)
                if n_nitrogen >= 2:
                    nucleobase_found = True
                    break
            if nucleobase_found:
                break
        if nucleobase_found:
            break
    if not nucleobase_found:
        return False, "No nucleobase found attached to the sugar (no neighboring aromatic ring with ≥2 nitrogen atoms)"
    
    # 4. Verify that a phosphate group is connected to the sugar ring.
    # We now attempt to find an exocyclic oxygen: a neighbor of a sugar ring atom that is not in the ring.
    # Then check if that oxygen is bonded to a phosphorus atom.
    phosphate_attached = False
    for atom_idx in sugar_ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in sugar_set:
                continue  # skip atoms inside the ring
            if nbr.GetAtomicNum() != 8:
                continue  # hoping for an oxygen atom
            # Check if this exocyclic oxygen is bonded to any phosphorus
            for oxy_nbr in nbr.GetNeighbors():
                if oxy_nbr.GetAtomicNum() == 15:
                    phosphate_attached = True
                    break
            if phosphate_attached:
                break
        if phosphate_attached:
            break
    
    # Also check if a phosphorus is directly bonded to a sugar atom (in case it is not separated by oxygen)
    if not phosphate_attached:
        for atom_idx in sugar_ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 15:
                    phosphate_attached = True
                    break
            if phosphate_attached:
                break

    if not phosphate_attached:
        return False, "No phosphate group attached to the sugar ring (expected a linkage via an oxygen substituent or directly)"
    
    return True, "Molecule contains a suitable sugar ring with an attached nucleobase and a phosphate group (via an exocyclic oxygen or direct bond) indicative of a nucleoside phosphate"

# Example usage when running as a script:
if __name__ == "__main__":
    test_smiles = "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N"  # 2'-deoxy-5-methyl-5'-cytidylic acid
    result, reason = is_nucleoside_phosphate(test_smiles)
    print(result, reason)