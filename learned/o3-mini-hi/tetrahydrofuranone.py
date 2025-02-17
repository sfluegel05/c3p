"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: Tetrahydrofuranone
Definition: Any oxolane (i.e. a 5‐membered ring with exactly one oxygen and four carbons)
that has at least one carbon atom in that ring bearing an exocyclic oxo substituent
(a C=O group, where the oxygen is not part of any ring).
Extra criteria to reduce false positives:
  • Each ring atom must not be aromatic.
  • Each ring atom must belong to no more than 2 rings.
  • The exocyclic oxygen must be terminal (degree = 1) and the ring carbon bearing it should have
    a total degree of 3 (two ring bonds plus the double bond to O).
"""

from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    
    A tetrahydrofuranone is defined as any molecule containing a 5‐membered ring (oxolane)
    having exactly one oxygen (and four carbons) in the ring, with at least one of the carbon atoms 
    in that ring carrying an exocyclic oxo substituent (a double‐bonded C=O group), where:
      - The ring is not aromatic.
      - No atom of the ring belongs to more than 2 rings.
      - The exocyclic oxygen is terminal (degree = 1) and not part of any ring.
      - The ring carbon bearing the C=O (the lactone carbonyl) shows the typical connectivity (degree equal to 3).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as tetrahydrofuranone, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Iterate over all rings present in the molecule.
    for ring in rings:
        # Only consider rings of size 5.
        if len(ring) != 5:
            continue
        
        # Check that the ring contains exactly 1 oxygen and 4 carbons and no other elements.
        oxygen_count = 0
        carbon_indices = []
        valid_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Reject if the atom is aromatic.
            if atom.GetIsAromatic():
                valid_ring = False
                break
            # Reject if atom belongs to >2 rings (overly fused system).
            if ring_info.NumAtomRings(idx) > 2:
                valid_ring = False
                break
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 8:
                oxygen_count += 1
            elif atomic_num == 6:
                carbon_indices.append(idx)
            else:
                valid_ring = False
                break
        if not valid_ring:
            continue
        # Must be exactly 1 oxygen and 4 carbons.
        if oxygen_count != 1 or len(carbon_indices) != 4:
            continue
        
        # Now, among the ring carbons, look for one that has an exocyclic C=O double bond.
        # We also require that the ring-carbon (bearing the carbonyl) has degree == 3.
        for c_idx in carbon_indices:
            c_atom = mol.GetAtomWithIdx(c_idx)
            # For a typical lactone, the carbonyl carbon should have exactly 3 bonds:
            # two connecting into the ring and the carbonyl bond.
            if c_atom.GetDegree() != 3:
                continue
            # Examine bonds on this ring carbon.
            for bond in c_atom.GetBonds():
                # We need a double bond to oxygen.
                if bond.GetBondType() != Chem.rdchem.BondType.DOUBLE:
                    continue
                neighbor = bond.GetOtherAtom(c_atom)
                # Candidate neighbor must be oxygen.
                if neighbor.GetAtomicNum() != 8:
                    continue
                # It must NOT be part of the candidate ring.
                if neighbor.GetIdx() in ring:
                    continue
                # The exocyclic oxygen must be terminal.
                if neighbor.GetDegree() != 1:
                    continue
                # When all criteria are met we accept the ring.
                return True, (
                    f"Found tetrahydrofuranone: ring (atom indices {ring}) has an oxo group "
                    f"attached to ring carbon index {c_idx} (exocyclic oxygen idx {neighbor.GetIdx()})."
                )
    # If no candidate ring fulfills all of the criteria.
    return False, "No tetrahydrofuranone substructure found."


# For testing purposes only:
if __name__ == "__main__":
    test_examples = [
        # True positives:
        ("N-butyryl-L-homoserine lactone", "CCCC(=O)N[C@@H]1CCOC1=O"),
        ("1,2,3,6-tetrahydrophthalic anhydride", "O=C1OC(=O)C2CC=CCC12"),
        ("3-[2-(phenylsulfonyloxy)ethyl]-5,5-dimethyldihydro-2(3H)-furanone", "CC1(C)CC(CCOS(=O)(=O)c2ccccc2)C(=O)O1"),
        ("(-)-5'-desmethylyatein", "C([C@@H]1[C@H](COC1=O)CC2=CC=C3OCOC3=C2)C4=CC(=C(C(=C4)O)OC)OC"),
        ("gamma-decalactone", "C1(OC(=O)CC1)CCCCCC"),
        ("L-homoserine lactone", "N[C@H]1CCOC1=O"),
        ("3-chloromethyl-5,5-dimethylbutyrolactone", "CC1(C)CC(CCl)C(=O)O1"),
        ("gamma-caprolactone", "CCC1CCC(=O)O1"),
        # Some examples known to be false positives in the previous attempt:
        ("Divaricoside (false positive test)", "CO[C@H]1C[C@H](O[C@H]2CC[C@@]3(C)[C@H](CC[C@@H]4[C@@H]3[C@H](O)C[C@]3(C)[C@H](CC[C@]43O)C)[H]")
    ]
    for name, smi in test_examples:
        res, reason = is_tetrahydrofuranone(smi)
        print(f"{name}:\n  Classified = {res}\n  Reason: {reason}\n")