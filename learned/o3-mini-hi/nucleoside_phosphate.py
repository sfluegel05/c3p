"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
#!/usr/bin/env python
"""
Classifies: Nucleoside phosphate
A nucleobase‐containing molecular entity that is a nucleoside in which 
one or more of the sugar hydroxy groups has been converted into a mono‐ or poly‐phosphate.
This includes both nucleotides and “non‐nucleotide” nucleoside phosphates.

Strategy:
  1. Parse the SMILES and require at least one phosphorus atom.
  2. Look for a candidate sugar ring: we require a non‐aromatic 5–membered ring with exactly 1 oxygen and 4 carbons.
  3. Check that at least one substituent off that sugar ring belongs to an aromatic 5– or 6–membered ring with ≥2 nitrogens (nucleobase).
  4. Verify that a phosphate is “attached” to the sugar by checking that either:
       – a non‐ring oxygen substituent of a sugar ring atom is directly bound to phosphorus, or
       – an exocyclic carbon (e.g. the 5′–CH₂ of ribose) has an oxygen substituent which in turn is bound to phosphorus.
       
If all conditions are met we return True.
Note: This is one possible implementation; there are many edge cases.
"""
from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple of (True, reason) if found, else (False, explanation)
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
    # 2. Look for a candidate sugar ring.
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Discard aromatic rings.
        if any(atom.GetIsAromatic() for atom in atoms):
            continue
        # Count oxygen and carbon atoms in ring.
        num_oxygens = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        num_carbons = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if num_oxygens == 1 and num_carbons == 4:
            sugar_ring = ring
            break
    if sugar_ring is None:
        return False, "No suitable sugar ring (non-aromatic 5‐membered ring with 1 O and 4 C) found"
    
    sugar_set = set(sugar_ring)
    
    # 3. Check for a nucleobase: at least one non‐ring neighbor of a sugar ring atom that is part of an aromatic ring (size 5 or 6) with ≥2 nitrogens.
    nucleobase_found = False
    for idx in sugar_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in sugar_set:
                continue
            # For each ring that this neighbor is in:
            for ring in ring_info.AtomRings():
                if nbr.GetIdx() not in ring:
                    continue
                if len(ring) not in (5, 6):
                    continue
                ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                # Require that all atoms in this ring are aromatic
                if not all(a.GetIsAromatic() for a in ring_atoms):
                    continue
                # Count nitrogen atoms
                n_nitrogen = sum(1 for a in ring_atoms if a.GetAtomicNum() == 7)
                if n_nitrogen >= 2:
                    nucleobase_found = True
                    break
            if nucleobase_found:
                break
        if nucleobase_found:
            break
    if not nucleobase_found:
        return False, "No nucleobase found attached to the sugar (no neighboring aromatic ring with ≥2 N atoms)"
    
    # 4. Check for a phosphate group attached to the sugar.
    phosphate_attached = False
    # First, for each sugar ring atom, check its immediate substituents:
    for idx in sugar_ring:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in sugar_set:
                continue  # Skip atoms in the ring.
            # Case 1: The neighbor is an oxygen that is directly bonded to phosphorus.
            if nbr.GetAtomicNum() == 8:
                for oxy_nbr in nbr.GetNeighbors():
                    if oxy_nbr.GetAtomicNum() == 15:
                        phosphate_attached = True
                        break
                if phosphate_attached:
                    break
            # Case 2: The neighbor is a carbon (e.g. exocyclic CH2) that itself might carry an oxygen bound to phosphorus.
            if nbr.GetAtomicNum() == 6:
                for sub in nbr.GetNeighbors():
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
            # Also handle the possibility that a sugar ring atom is directly bound to phosphorus.
            if nbr.GetAtomicNum() == 15:
                phosphate_attached = True
                break
        if phosphate_attached:
            break
    if not phosphate_attached:
        return False, "No phosphate group attached to the sugar ring (expected via an intervening oxygen or direct bond)"
    
    return True, "Molecule contains a suitable furanose sugar with an attached nucleobase and a phosphate (directly or via a CH2 substituent) indicative of a nucleoside phosphate"

# Example usage when running as a script:
if __name__ == "__main__":
    # Test example – one of the provided molecules.
    test_smiles = "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N"  # 2'-deoxy-5-methyl-5'-cytidylic acid (false negative previously)
    result, reason = is_nucleoside_phosphate(test_smiles)
    print(result, reason)