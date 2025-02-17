"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: Tetrahydrofuranone
Definition: Any oxolane (i.e. a 5‐membered ring with exactly one oxygen and four carbons)
that has at least one carbon atom in that ring bearing an exocyclic oxo substituent (a C=O group).
We now add extra checks to reduce false positives by:
  • Rejecting candidate rings with any aromatic atoms.
  • Ensuring that the exocyclic oxygen (in the carbonyl) is not itself part of a ring and has degree 1.
  • Avoiding overly fused rings by ensuring that no atom of the candidate ring belongs to more than 2 rings.
"""

from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    
    A tetrahydrofuranone is defined as any molecule containing a 5‐membered ring (oxolane)
    having exactly one oxygen (and four carbons) in the ring, with at least one of the carbon atoms 
    in that ring carrying an exocyclic oxo substituent (a double‐bonded C=O group, where the O is not in the ring).
    
    Extra criteria:
      - The ring atoms must not be aromatic.
      - The exocyclic oxygen must be terminal (degree 1).
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
    
    # Iterate over all rings in the molecule.
    for ring in rings:
        if len(ring) != 5:
            continue
        # Separate ring atoms by element type and collect carbon atom indices.
        oxygen_count = 0
        carbon_indices = []
        skip_ring = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Reject candidate if any atom is aromatic.
            if atom.GetIsAromatic():
                skip_ring = True
                break
            # Also, avoid rings where atoms are overly fused (i.e. belong to >2 rings)
            if ring_info.NumAtomRings(idx) > 2:
                skip_ring = True
                break
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                carbon_indices.append(idx)
            else:
                # If there is an atom that is not carbon or oxygen, this ring is not an oxolane.
                skip_ring = True
                break
        if skip_ring:
            continue
        if oxygen_count != 1 or len(carbon_indices) != 4:
            continue
        
        # For each ring carbon, check if it bears an exocyclic oxo substituent.
        for c_idx in carbon_indices:
            c_atom = mol.GetAtomWithIdx(c_idx)
            for bond in c_atom.GetBonds():
                # Look only at bonds that are double bonds.
                if bond.GetBondTypeAsDouble() != 2.0:
                    continue
                # Identify the neighbor on the other side.
                neighbor = bond.GetOtherAtom(c_atom)
                # We want an oxygen atom that is not part of the ring (exocyclic).
                if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in ring:
                    # It must be a terminal oxygen (degree 1) to avoid misinterpreting bridging oxygens.
                    if neighbor.GetDegree() == 1:
                        return True, f"Found tetrahydrofuranone: ring (atom indices {ring}) has an oxo group at ring carbon index {c_idx} (exocyclic oxygen idx {neighbor.GetIdx()})."
        # End for each candidate ring carbon.
    # End for each ring.
    
    return False, "No tetrahydrofuranone substructure found."

# Example usage (for testing purposes only):
if __name__ == "__main__":
    examples = [
        # True positives:
        ("N-butyryl-L-homoserine lactone", "CCCC(=O)N[C@@H]1CCOC1=O"),
        ("1,2,3,6-tetrahydrophthalic anhydride", "O=C1OC(=O)C2CC=CCC12"),
        ("3-[2-(phenylsulfonyloxy)ethyl]-5,5-dimethyldihydro-2(3H)-furanone", "CC1(C)CC(CCOS(=O)(=O)c2ccccc2)C(=O)O1"),
        ("(-)-5'-desmethylyatein", "C([C@@H]1[C@H](COC1=O)CC2=CC=C3OCOC3=C2)C4=CC(=C(C(=C4)O)OC)OC"),
        ("gamma-decalactone", "C1(OC(=O)CC1)CCCCCC"),
        ("L-homoserine lactone", "N[C@H]1CCOC1=O"),
        ("gamma-caprolactone", "CCC1CCC(=O)O1"),
        # A false positive test example (should return False):
        ("Divaricoside (false positive test)", "CO[C@H]1C[C@H](O[C@H]2CC[C@@]3(C)[C@H](CC[C@@H]4[C@@H]3[C@H](O)C[C@]3(C)[C@H](CC[C@]43O)C)[C@@]1([H])")
    ]
    for name, smi in examples:
        res, reason = is_tetrahydrofuranone(smi)
        print(f"{name}:\n  Classified = {res}\n  Reason: {reason}\n")