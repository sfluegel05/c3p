"""
Classifies: CHEBI:28863 flavanones
"""
#!/usr/bin/env python3
"""
Classifies: flavanones – Defined as members of the class of flavans with a 
3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.
This implementation performs ring-based analysis to require:
  • a six-membered non-aromatic ring that has exactly one oxygen,
  • a carbonyl group (a C=O where the oxygen is exocyclic),
  • fusion (sharing exactly two atoms) with an aromatic six‐membered ring,
  • and at least one atom (other than the carbonyl carbon) of the candidate ring bears
    an external aromatic substituent (proxy for the “2-aryl” group).
"""

from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.

    The algorithm works as follows:
      1. Parses the molecule and extracts all ring systems.
      2. For each six-membered ring that is non-aromatic and contains exactly one oxygen,
         it checks for a carbonyl group (a double-bonded oxygen attached exocyclically).
      3. It then looks for a fused aromatic ring (a benzene ring) that shares exactly two atoms
         with the candidate six‐membered ring.
      4. Finally, it requires that one of the non–carbonyl atoms in the candidate ring bears an 
         external aromatic substituent (indicative of the "2-aryl" group).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): True and a reason if classified as a flavanone, else False with a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Pre-collect all aromatic rings of size 6 (candidate benzene rings)
    aromatic_6_rings = []
    for ring in ring_info:
        if len(ring) == 6:
            # Check if all atoms in this ring are aromatic
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                aromatic_6_rings.append(set(ring))
    
    # Loop over all candidate 6-membered rings in the molecule
    for candidate in ring_info:
        if len(candidate) != 6:
            continue  # focus only on six-membered rings
        
        # Get the atoms in this candidate ring.
        candidate_atoms = [mol.GetAtomWithIdx(idx) for idx in candidate]
        # The dihydropyranone ring should be non-aromatic 
        if any(atom.GetIsAromatic() for atom in candidate_atoms):
            continue
        
        # Count oxygen atoms (should be exactly one in the ring)
        oxygens_in_ring = [atom for atom in candidate_atoms if atom.GetAtomicNum() == 8]
        if len(oxygens_in_ring) != 1:
            continue
        
        # Try to find a carbonyl carbon in the ring: a carbon (atomic num 6) 
        # that has at least one double bond to an oxygen that is not in the candidate ring.
        carbonyl_atom_idx = None
        for atom in candidate_atoms:
            if atom.GetAtomicNum() != 6:
                continue
            for bond in atom.GetBonds():
                # Check if bond is a double bond
                if bond.GetBondTypeAsDouble() == 2:
                    # Identify the neighbor atom
                    nbr = bond.GetOtherAtom(atom)
                    # It should be oxygen and not part of the candidate ring.
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in candidate:
                        carbonyl_atom_idx = atom.GetIdx()
                        break
            if carbonyl_atom_idx is not None:
                break
        if carbonyl_atom_idx is None:
            continue  # no carbonyl group found in this ring
        
        # Next, check that the candidate ring is fused to a benzene ring.
        fused_aromatic_ring = None
        for arom_ring in aromatic_6_rings:
            # Fused rings share at least 2 atoms. 
            if len(arom_ring.intersection(set(candidate))) == 2:
                fused_aromatic_ring = arom_ring
                break
        if fused_aromatic_ring is None:
            continue  # no fused benzene ring found
        
        # Finally, check that one of the non-carbonyl atoms in candidate ring
        # bears an external aromatic substituent (that is not in the fused aromatic ring).
        aryl_substituent_found = False
        for idx in candidate:
            if idx == carbonyl_atom_idx:
                continue  # skip the carbonyl carbon
            atom = mol.GetAtomWithIdx(idx)
            # Look at neighbors not in the candidate ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in candidate:
                    continue
                # We want the substituent neighbor to be aromatic, in a ring, 
                # and not simply part of the fused benzene ring.
                if nbr.GetIsAromatic() and nbr.IsInRing():
                    # If the neighbor belongs to a ring that is not exactly the fused aromatic ring...
                    neighbor_in_fused = False
                    for ring in ring_info:
                        if nbr.GetIdx() in ring and set(ring) == fused_aromatic_ring:
                            neighbor_in_fused = True
                            break
                    if not neighbor_in_fused:
                        aryl_substituent_found = True
                        break
            if aryl_substituent_found:
                break
        if not aryl_substituent_found:
            continue  # no external aryl substituent detected on candidate ring
        
        # If all criteria are met, we classify the molecule as flavanone.
        return True, ("Molecule contains the flavanone core: a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton "
                      "with a dihydropyranone fused to a benzene ring and an external aryl substituent.")
    
    # If no candidate ring qualifies, return false.
    return False, "Flavanone core structure (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one) not found"

# Example usage for testing (if executed as main)
if __name__ == "__main__":
    test_examples = {
        "eriodictyol": "Oc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1",
        "sakuranetin": "COc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)cc1",
        "(R)-naringenin": "Oc1ccc(cc1)[C@H]1CC(=O)c2c(O)cc(O)cc2O1"
    }
    for name, smi in test_examples.items():
        res, reason = is_flavanones(smi)
        print(f"{name}: {res} ({reason})")