"""
Classifies: CHEBI:28863 flavanones
"""
#!/usr/bin/env python3
"""
Classifies: flavanones – Defined as members of the class of flavans with a 
3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.
This implementation does not rely solely on a single SMARTS pattern; instead,
it uses ring analysis to detect:
  • a six-membered, non-aromatic heterocycle containing exactly one oxygen atom,
  • a carbonyl group (C(=O)) in that ring, and
  • at least one ring atom (a candidate for C-2) that bears a substituent not in the ring which is aromatic.
Additionally, it checks that the heterocycle is fused to an aromatic ring.
"""

from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.

    The algorithm searches for a six-membered, non-aromatic ring (the dihydropyran part)
    that (a) contains exactly one oxygen; (b) contains a carbonyl carbon (double-bonded O);
    (c) is fused to an aromatic ring (the benzene ring of the benzopyran portion); and
    (d) has at least one non–carbonyl ring carbon that bears an external aromatic substituent
        (as a proxy for the "2-aryl" group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a flavanone, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo().AtomRings()
    # For each ring, check if it could be the heterocyclic (dihydro-pyranone) ring
    for ring in ring_info:
        if len(ring) != 6:
            continue  # we focus only on six-membered rings
        
        # Get the atoms in this ring
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check that none (or very few) of the atoms are aromatic; in a dihydro ring they should be non-aromatic.
        if any(atom.GetIsAromatic() for atom in ring_atoms):
            # Skip rings that are fully aromatic.
            continue
        
        # Count oxygen atoms in the ring. For flavanone, there should be exactly one (O at position 1).
        oxygens_in_ring = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
        if len(oxygens_in_ring) != 1:
            continue

        # Look for a carbonyl group within the ring.
        carbonyl_atom = None
        for atom in ring_atoms:
            if atom.GetAtomicNum() != 6:
                continue
            # A carbonyl carbon should have at least one double bond to oxygen.
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() == 2:  # double bond
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8:
                        # Found a carbonyl. We record this atom.
                        carbonyl_atom = atom
                        break
            if carbonyl_atom is not None:
                break
        if carbonyl_atom is None:
            continue  # no carbonyl found in this ring
        
        # Check that the candidate ring is fused to an aromatic ring.
        # Look at each atom in the ring: if any neighbor outside the ring is aromatic and belongs to a small (likely benzene) ring, we take that as evidence.
        fused_aromatic = False
        for atom in ring_atoms:
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetIsAromatic() and nbr.IsInRing():
                    fused_aromatic = True
                    break
            if fused_aromatic:
                break
        if not fused_aromatic:
            continue  # candidate heterocycle not fused to an aromatic ring
        
        # Check that one of the non-carbonyl carbons in the ring bears an external aromatic substituent.
        # We call this candidate the "C-2" position.
        aryl_substituent_found = False
        for atom in ring_atoms:
            # Consider only carbons that are not the carbonyl carbon
            if atom.GetAtomicNum() != 6 or atom.GetIdx() == carbonyl_atom.GetIdx():
                continue
            # Look for neighbors outside of the ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # If the substituent neighbor is aromatic and is in a ring, we take it as evidence.
                if nbr.GetIsAromatic() and nbr.IsInRing():
                    aryl_substituent_found = True
                    break
            if aryl_substituent_found:
                break
        if not aryl_substituent_found:
            continue  # no external aryl substituent detected on a candidate ring carbon
        
        # If we reached here, the candidate ring meets all our criteria.
        return True, "Molecule contains the flavanone core: a benzopyran-4-one system with an aryl substituent at C-2"

    # If no candidate ring qualifies, return false with reason.
    return False, "Flavanone core structure (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one) not found"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Examples taken from the supplied list (few examples shown for brevity)
    examples = {
        "eriodictyol": "Oc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1",
        "sakuranetin": "COc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)cc1",
        "(R)-naringenin": "Oc1ccc(cc1)[C@H]1CC(=O)c2c(O)cc(O)cc2O1",
        "hesperetin": "COc1ccc(cc1O)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1",
        "diosmin": "COc1ccc(cc1O)-c1cc(=O)c2c(O)cc(O)c2o1"
    }
    
    for name, smi in examples.items():
        result, reason = is_flavanones(smi)
        print(f"{name}: {result} ({reason})")