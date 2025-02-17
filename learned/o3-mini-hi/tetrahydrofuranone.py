"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: Tetrahydrofuranone
Definition: Any oxolane (i.e. a 5‐membered ring with exactly one oxygen and four carbons)
that has at least one carbon atom in that ring bearing an exocyclic oxo substituent (a C=O group).
Extra criteria have been imposed to reduce false positives:
  • Reject candidate rings with any aromatic atoms.
  • Ensure that the exocyclic oxygen (in the carbonyl) is not itself part of any ring and is terminal (degree = 1).
  • Avoid overly fused rings by ensuring that no atom of the candidate ring belongs to more than 2 rings.
"""

from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    
    A tetrahydrofuranone is defined as any molecule containing a 5‐membered ring (oxolane)
    having exactly one oxygen (and four carbons) in the ring, with at least one of the carbon atoms 
    in that ring carrying an exocyclic oxo substituent (a double‐bonded C=O group, where the O is not part of the ring).
    
    Extra criteria:
      - The ring atoms must not be aromatic.
      - The exocyclic oxygen must be terminal (degree = 1).
      - To avoid candidate rings that are part of a highly fused system, no atom in the candidate
        ring should participate in more than 2 rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as tetrahydrofuranone, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Iterate through each ring of size 5
    for ring in rings:
        if len(ring) != 5:
            continue
        
        oxygen_count = 0
        carbon_indices = []
        valid_ring = True  # assume it is a candidate unless it fails one of the tests
        
        # Check each atom in the ring
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Reject if atom is aromatic
            if atom.GetIsAromatic():
                valid_ring = False
                break
            # Reject overly fused atoms (an atom present in >2 rings)
            if ring_info.NumAtomRings(idx) > 2:
                valid_ring = False
                break
            # Only consider carbon and oxygen atoms
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                carbon_indices.append(idx)
            else:
                valid_ring = False
                break
        if not valid_ring:
            continue
        
        # The candidate ring must have exactly 1 oxygen and 4 carbons.
        if oxygen_count != 1 or len(carbon_indices) != 4:
            continue
        
        # For each carbon in the candidate ring, check its bonds for an exocyclic C=O double bond.
        for c_idx in carbon_indices:
            c_atom = mol.GetAtomWithIdx(c_idx)
            for bond in c_atom.GetBonds():
                # Only look at bonds with double bond type (use bond type comparison)
                if bond.GetBondType() != Chem.rdchem.BondType.DOUBLE:
                    continue
                neighbor = bond.GetOtherAtom(c_atom)
                # The neighbor must be oxygen and not be a member of this ring.
                if neighbor.GetAtomicNum() != 8 or neighbor.GetIdx() in ring:
                    continue
                # The exocyclic oxygen should be terminal (degree exactly =1).
                if neighbor.GetDegree() != 1:
                    continue

                # Found a candidate exocyclic carbonyl on a ring carbon.
                return True, (
                    f"Found tetrahydrofuranone: ring (atom indices {ring}) has an oxo group "
                    f"at ring carbon index {c_idx} (exocyclic oxygen idx {neighbor.GetIdx()})."
                )
    # End candidate ring search.
    return False, "No tetrahydrofuranone substructure found."

# For testing purposes only:
if __name__ == "__main__":
    examples = [
        # True positives:
        ("N-butyryl-L-homoserine lactone", "CCCC(=O)N[C@@H]1CCOC1=O"),
        ("1,2,3,6-tetrahydrophthalic anhydride", "O=C1OC(=O)C2CC=CCC12"),
        ("3-[2-(phenylsulfonyloxy)ethyl]-5,5-dimethyldihydro-2(3H)-furanone", "CC1(C)CC(CCOS(=O)(=O)c2ccccc2)C(=O)O1"),
        ("(-)-5'-desmethylyatein", "C([C@@H]1[C@H](COC1=O)CC2=CC=C3OCOC3=C2)C4=CC(=C(C(=C4)O)OC)OC"),
        ("gamma-decalactone", "C1(OC(=O)CC1)CCCCCC"),
        ("L-homoserine lactone", "N[C@H]1CCOC1=O"),
        ("3-chloromethyl-5,5-dimethylbutyrolactone", "CC1(C)CC(CCl)C(=O)O1"),
        ("gamma-caprolactone", "CCC1CCC(=O)O1"),
        # One false positive candidate (should return False):
        ("Divaricoside (false positive test)", "CO[C@H]1C[C@H](O[C@H]2CC[C@@]3(C)[C@H](CC[C@@H]4[C@@H]3[C@H](O)C[C@]3(C)[C@H](CC[C@]43O)C)[H]")
    ]
    for name, smi in examples:
        res, reason = is_tetrahydrofuranone(smi)
        print(f"{name}:\n  Classified = {res}\n  Reason: {reason}\n")