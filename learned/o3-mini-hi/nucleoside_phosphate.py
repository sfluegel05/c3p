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
  2. Check for the existence of at least one phosphorus atom.
  3. Identify a candidate sugar ring: a non‐aromatic 5‐membered ring with exactly one oxygen and four carbons.
  4. Verify that one atom in that sugar ring bonds to a nucleobase – here defined as a neighboring ring (5- or 6-membered) that is aromatic and contains at least two nitrogen atoms.
  5. Verify that a phosphate group is attached to the sugar ring by checking that a phosphorus atom bonds to an oxygen that belongs to the sugar ring.
  
If all tests pass, classify the molecule as a nucleoside phosphate.
"""
from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True and a success reason if classification is successful,
                     False with an explanation otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Basic check: Should contain at least one phosphorus atom.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_atoms:
        return False, "No phosphorus atoms found; unlikely to be a nucleoside phosphate"
    
    # 2. Identify a candidate sugar ring:
    ring_info = mol.GetRingInfo()
    sugar_ring = None
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # For sugar rings we expect them to be non-aromatic.
        if any(atom.GetIsAromatic() for atom in atoms):
            continue
        num_oxygens = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        num_carbons = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if num_oxygens == 1 and num_carbons == 4:
            sugar_ring = ring
            break
    if sugar_ring is None:
        return False, "No suitable sugar ring (non-aromatic 5-membered ring with 1 oxygen and 4 carbons) found"
    
    # 3. Identify a nucleobase attached to the sugar.
    # We look for a neighbor of any sugar atom that is part of an aromatic ring (5- or 6-membered)
    # with at least 2 nitrogen atoms. This ensures that the nucleobase is linked to the sugar.
    nucleobase_found = False
    sugar_set = set(sugar_ring)
    for idx in sugar_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in sugar_set:
                continue  # already in sugar; skip
            # Now, check if this neighbor participates in a ring that could be a nucleobase.
            for ring in ring_info.AtomRings():
                # We require that this neighbor is in a ring that is 5 or 6 atoms in size.
                if nbr_idx not in ring or len(ring) not in (5,6):
                    continue
                ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                # Check that the ring is aromatic.
                if not all(a.GetIsAromatic() for a in ring_atoms):
                    continue
                # Count nitrogen atoms.
                n_nitrogen = sum(1 for a in ring_atoms if a.GetAtomicNum() == 7)
                if n_nitrogen >= 2:
                    nucleobase_found = True
                    break
            if nucleobase_found:
                break
        if nucleobase_found:
            break
    if not nucleobase_found:
        return False, "No nucleobase found attached to the sugar ring (aromatic ring with ≥2 nitrogen atoms)"
    
    # 4. Verify that a phosphate group is connected to the sugar ring.
    # Rather than using a fixed SMARTS, we look at phosphorus atoms and see if any of their oxygen neighbors are in the sugar.
    phosphate_attached = False
    for pAtom in phosphorus_atoms:
        for nbr in pAtom.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetIdx() in sugar_set:
                phosphate_attached = True
                break
        if phosphate_attached:
            break
    if not phosphate_attached:
        return False, "No phosphate group attached to the sugar ring (expected P-O linkage)"
    
    # If all tests pass, we classify the molecule as a nucleoside phosphate.
    return True, "Molecule contains a sugar ring with an attached nucleobase and phosphate group indicative of a nucleoside phosphate"

# Example usage:
if __name__ == "__main__":
    # Example: 2'-deoxy-5-methyl-5'-cytidylic acid
    test_smiles = "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N"
    result, reason = is_nucleoside_phosphate(test_smiles)
    print(result, reason)