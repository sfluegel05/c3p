"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: Mucopolysaccharide (glycosaminoglycan)
Definition: “Any of the group of polysaccharides composed of alternating units from uronic acids
and glycosamines, and commonly partially esterified with sulfuric acid.”
This improved implementation uses a heuristic that:
  1. Parses the SMILES into a molecule.
  2. Loops over six‐membered rings that have only carbon and oxygen atoms (typical pyranose rings).
  3. For each candidate ring, examines exocyclic substituents:
       • A directly attached nitrogen (ignoring amide-like functionality) flags the ring as glycosamine–like.
       • A substituent carbon that bears a double‐bonded oxygen plus at least one more oxygen (putative carboxyl group) flags it as uronic acid–like.
       • If both features appear the ring is skipped as ambiguous.
  4. It then searches for glycosidic connectivity (an exocyclic oxygen linking atoms from two candidate rings).
  5. Finally, the code requires a minimum total number of candidate rings,
     and at least two of each unambiguous type (with nearly balanced counts), as well as a minimum molecular weight.
Note: This is a heuristic approach; many complex carbohydrate motifs may be missed or spuriously flagged.
If one judges the challenge too hard to reliably classify, returning (None, None) is acceptable.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide (glycosaminoglycan) based on a sugar‐ring heuristic.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a mucopolysaccharide, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Gather candidate rings. We look only for six-membered rings (pyranose-like)
    rings = mol.GetRingInfo().AtomRings()
    candidate_rings = []  # Each element will be (ring, ring_type) where ring_type is either "glyco" or "uronic"
    
    # Loop over rings and filter for those with acceptable composition:
    for ring in rings:
        # Only consider rings of size six.
        if len(ring) != 6:
            continue
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Require that the ring atoms are only carbon (6) or oxygen (8)
        if any(atom.GetAtomicNum() not in (6, 8) for atom in ring_atoms):
            continue
        # Require exactly one oxygen atom in the ring (typical for pyranose)
        if sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8) != 1:
            continue

        # Initialize flags for this candidate ring:
        glyco_flag = False
        uronic_flag = False

        # Look at each ring carbon’s exocyclic substituents
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # We consider only carbons in the ring for substituent analysis.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # Skip atoms that are in the ring.
                # Check for glycosamine indicator: a directly attached nitrogen.
                if nbr.GetAtomicNum() == 7:
                    # We try to avoid flagging an amide nitrogen by checking that the neighbor
                    # is not directly connected to a carbonyl (this is a simple check).
                    has_carbonyl = False
                    for nn in nbr.GetNeighbors():
                        if nn.GetAtomicNum() == 6:
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                                has_carbonyl = True
                    if not has_carbonyl:
                        glyco_flag = True

                # Check for uronic acid indicator: a substituent carbon that carries a carboxyl group.
                if nbr.GetAtomicNum() == 6:
                    ext_neighbors = [n for n in nbr.GetNeighbors() if n.GetIdx() not in (ring + [atom.GetIdx()])]
                    ox_count = 0
                    dbl_ox = 0
                    for ext in ext_neighbors:
                        if ext.GetAtomicNum() == 8:
                            ox_count += 1
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), ext.GetIdx())
                            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                                dbl_ox += 1
                    # A minimal carboxyl group: at least one double bond and at least one extra oxygen.
                    if dbl_ox >= 1 and ox_count >= 2:
                        uronic_flag = True
        
        # If a ring shows both indicators, skip as ambiguous.
        if glyco_flag and uronic_flag:
            continue
        # Otherwise add the ring if it is flagged as one or the other.
        if glyco_flag or uronic_flag:
            ring_type = "glyco" if glyco_flag else "uronic"
            candidate_rings.append((ring, ring_type))
    
    # Count rings of each type.
    glyco_count = sum(1 for ring, typ in candidate_rings if typ == "glyco")
    uronic_count = sum(1 for ring, typ in candidate_rings if typ == "uronic")
    total_rings = len(candidate_rings)
    
    # Check for glycosidic connectivity: does an exocyclic oxygen (from a ring carbon) connect two candidate rings?
    connectivity = False
    for i in range(total_rings):
        ring_i = candidate_rings[i][0]
        for j in range(i+1, total_rings):
            ring_j = candidate_rings[j][0]
            # Loop over each carbon atom in ring_i and its exocyclic oxygens.
            for idx in ring_i:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring_i:
                        for nbr2 in nbr.GetNeighbors():
                            if nbr2.GetIdx() in ring_j:
                                connectivity = True
                                break
                    if connectivity:
                        break
                if connectivity:
                    break
            if connectivity:
                break
        if connectivity:
            break
    
    summary = (f"Detected {uronic_count} uronic acid-like and {glyco_count} glycosamine-like rings "
               f"across {total_rings} candidate pyranose ring(s).")
    summary += " " + ("Glycosidic connectivity detected." if connectivity else "No glycosidic connectivity detected.")
    
    # Apply heuristic thresholds:
    if total_rings < 4:
        return False, f"Not enough sugar-like rings detected (found {total_rings}, need at least 4)."
    if uronic_count < 2 or glyco_count < 2:
        return False, (f"Repeating unit counts too low: {uronic_count} uronic acid-like and "
                       f"{glyco_count} glycosamine-like rings (need at least 2 of each).")
    if abs(uronic_count - glyco_count) > 1:
        return False, (f"Repeating unit counts imbalanced: {uronic_count} uronic acid-like vs. "
                       f"{glyco_count} glycosamine-like rings (difference > 1).")
    
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 1000:
        return False, f"Molecular weight too low for a typical mucopolysaccharide ({mw:.1f} Da)."
    
    return True, summary

# Example usage:
if __name__ == "__main__":
    # Replace the following SMILES with one of your test examples.
    example_smiles = "O=C1NC(C(=O)NCCCC(=O)O)CC1"  
    result, reason = is_mucopolysaccharide(example_smiles)
    print("Result:", result)
    print("Reason:", reason)