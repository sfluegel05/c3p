"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: Mucopolysaccharide (glycosaminoglycan)
Definition: “Any of the group of polysaccharides composed of alternating units from uronic acids
and glycosamines, and commonly partially esterified with sulfuric acid.”
This implementation uses a heuristic that:
  1. Parses the SMILES into a molecule.
  2. Loops over six‐membered rings that have only carbon and oxygen atoms (typical pyranose rings).
  3. For each candidate ring, examines exocyclic substituents:
       • A directly attached nitrogen (but ignoring typical amide nitrogens) flags the ring as glycosamine–like.
       • A substituent carbon that carries a double‐bonded oxygen plus at least one extra oxygen (putative carboxyl)
         flags it as uronic acid–like.
       • Rings that show both features are skipped as ambiguous.
  4. It then searches for glycosidic connectivity: an exocyclic oxygen linking two candidate rings.
  5. Finally, the code requires a minimum total number of candidate rings, at least two rings of each unambiguous type,
     and a minimum molecular weight.
Note: It is a heuristic and many complex carbohydrates may be missed or spuriously flagged.
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
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Gather candidate rings. We focus on six-membered rings (pyranose-like).
    rings = mol.GetRingInfo().AtomRings()
    candidate_rings = []  # Each item: (ring (tuple of atom indices), ring_type [either "glyco" or "uronic"])
    
    # Process each ring.
    for ring in rings:
        # Only consider rings of size six.
        if len(ring) != 6:
            continue
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Ensure the ring is composed only of carbon (atomic 6) and oxygen (atomic 8).
        if any(atom.GetAtomicNum() not in (6, 8) for atom in ring_atoms):
            continue
        # Require exactly one oxygen atom in the ring (typical for pyranose rings).
        if sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8) != 1:
            continue

        # Flags for this candidate ring.
        glyco_flag = False
        uronic_flag = False

        # For each carbon in the ring, examine its exocyclic neighbors.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbon atoms.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # Skip atoms inside the ring.
                # Check for glycosamine indicator: directly attached nitrogen without a clear amide (carbonyl) connectivity.
                if nbr.GetAtomicNum() == 7:
                    has_carbonyl = False
                    for nn in nbr.GetNeighbors():
                        if nn.GetAtomicNum() == 6:
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                                has_carbonyl = True
                    if not has_carbonyl:
                        glyco_flag = True

                # Check for uronic acid indicator: look for a substituent carbon (neighbor) carrying a carboxyl group.
                if nbr.GetAtomicNum() == 6:
                    # Here we must ensure that we are not including atoms already in the ring.
                    # Since ring is a tuple and atom.GetIdx() is an int, we make a concatenated tuple.
                    exclude_atoms = tuple(ring) + (atom.GetIdx(),)
                    ext_neighbors = [n for n in nbr.GetNeighbors() if n.GetIdx() not in exclude_atoms]
                    ox_count = 0
                    dbl_ox = 0
                    for ext in ext_neighbors:
                        if ext.GetAtomicNum() == 8:
                            ox_count += 1
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), ext.GetIdx())
                            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                                dbl_ox += 1
                    # A minimal carboxyl group: at least one double-bonded oxygen and at least one extra oxygen.
                    if dbl_ox >= 1 and ox_count >= 2:
                        uronic_flag = True
        
        # Skip rings that show both indicators as ambiguous.
        if glyco_flag and uronic_flag:
            continue
        # Add the ring if it is flagged as one unambiguous type.
        if glyco_flag or uronic_flag:
            ring_type = "glyco" if glyco_flag else "uronic"
            candidate_rings.append((ring, ring_type))
    
    # Count rings of each type.
    glyco_count = sum(1 for ring, typ in candidate_rings if typ == "glyco")
    uronic_count = sum(1 for ring, typ in candidate_rings if typ == "uronic")
    total_rings = len(candidate_rings)
    
    # Check for glycosidic connectivity:
    # We examine if a directly attached oxygen from a carbon in one candidate ring connects to an atom in another ring.
    connectivity = False
    for i in range(total_rings):
        ring_i = candidate_rings[i][0]
        for j in range(i+1, total_rings):
            ring_j = candidate_rings[j][0]
            # Loop over each carbon atom in ring_i.
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
    
    # Apply heuristic thresholds.
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
    # For testing, one can replace this with any of the provided SMILES examples.
    example_smiles = "O=C1O[C@@H](C=C[C@@H](CC[C@]23O[C@H](C=4C(=C(NC(CC=C1C)=O)C=C(O)C4)O2)[C@@H](C)C(C3)=O)CC)[C@@H](O)C=C(C)C"
    result, reason = is_mucopolysaccharide(example_smiles)
    print("Result:", result)
    print("Reason:", reason)