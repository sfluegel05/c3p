"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: Mucopolysaccharide (glycosaminoglycan)
Definition: “Any of the group of polysaccharides composed of alternating units from uronic acids
and glycosamines, and commonly partially esterified with sulfuric acid.”
This implementation uses a heuristic that:
  1. Parses the SMILES into a molecule.
  2. Searches for candidate pyranose (6-membered) rings that contain exactly one oxygen (as in typical sugars)
     and all other ring atoms are carbons.
  3. For each candidate ring, examines substituents on its ring carbons that are not in the ring.
     - If a directly attached atom is nitrogen, the ring is marked as glycosamine–like.
     - If a substituent carbon bears at least one double‐bonded oxygen (and at least one more oxygen)
       it is taken as evidence for a carboxyl group (a uronic acid indicator).
     Rings showing both features are considered ambiguous and skipped.
  4. In addition, the code checks for glycosidic connectivity – i.e. if an exocyclic oxygen from one candidate
     ring is bound to an atom lying in another candidate ring.
  5. Finally, the code checks that the overall counts are high enough (at least 4 candidate rings, at least 2 of 
     each type with nearly balanced numbers) and that the overall molecular weight is high (arbitrarily set here).
  
Note: This algorithm is heuristic and may have both false positives and false negatives.
If one judges that this challenge is too difficult to reliably capture, one may return (None, None).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide (glycosaminoglycan) based on a sugar-ring heuristic.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if classified as mucopolysaccharide, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First pass: identify candidate pyranose rings.
    rings = mol.GetRingInfo().AtomRings()
    candidate_rings = []  # list of tuples: (ring, is_glyco, is_uronic)
    
    for ring in rings:
        # Only consider 6-membered rings (typical for pyranoses)
        if len(ring) != 6:
            continue
        
        # Get the atoms in this ring.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Require that every atom in the ring is either carbon (atomic#6) or oxygen (atomic#8)
        if any(atom.GetAtomicNum() not in (6, 8) for atom in ring_atoms):
            continue
        
        # Require that the ring has exactly one oxygen atom (like a typical pyranose ring)
        ring_oxy = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        if ring_oxy != 1:
            continue
        
        # Initialize flags (each ring is ideally either glycosamine–like or uronic acid–like)
        glyco_flag = False
        uronic_flag = False
        
        # For each ring carbon, check its exocyclic substituents.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Process only carbons (we expect the sugar carbons to bear revealing groups)
            if atom.GetAtomicNum() != 6:
                continue
            
            for nbr in atom.GetNeighbors():
                # Skip atoms that are in the ring.
                if nbr.GetIdx() in ring:
                    continue
                
                # Check for glycosamine indicator: directly attached nitrogen.
                if nbr.GetAtomicNum() == 7:
                    glyco_flag = True
                
                # Check for uronic acid indicator: a substituent carbon that bears a carboxyl group.
                if nbr.GetAtomicNum() == 6:
                    # Look at substituents of this exocyclic carbon (exclude the ring attachment)
                    ext_neighbors = [n for n in nbr.GetNeighbors() if n.GetIdx() not in ring and n.GetIdx() != atom.GetIdx()]
                    ox_count = 0
                    double_count = 0
                    for ext in ext_neighbors:
                        if ext.GetAtomicNum() == 8:
                            ox_count += 1
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), ext.GetIdx())
                            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                                double_count += 1
                    # Carboxyl group typical features: at least one double‐bonded oxygen and at least one additional oxygen.
                    if double_count >= 1 and ox_count >= 2:
                        uronic_flag = True
        
        # If a ring shows both glycosamine and uronic acid features simultaneously, skip it as ambiguous.
        if glyco_flag and uronic_flag:
            continue
        
        # Accept the ring only if it is clearly one type.
        if glyco_flag or uronic_flag:
            candidate_rings.append((ring, glyco_flag, uronic_flag))
    
    # Count the candidate rings by type.
    glyco_count = sum(1 for ring, is_glyco, is_uronic in candidate_rings if is_glyco)
    uronic_count = sum(1 for ring, is_glyco, is_uronic in candidate_rings if is_uronic)
    total_sugar_rings = len(candidate_rings)
    
    # Check connectivity: look for evidence of glycosidic bonds between candidate rings.
    connectivity = 0
    # Loop over each pair of candidate rings.
    for i in range(len(candidate_rings)):
        ring_i = candidate_rings[i][0]
        for j in range(i+1, len(candidate_rings)):
            ring_j = candidate_rings[j][0]
            # For each carbon in ring i, check its exocyclic oxygen substituents.
            found_connection = False
            for idx in ring_i:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue
                for nbr in atom.GetNeighbors():
                    # Only consider oxygen substituents that are not part of ring_i.
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring_i:
                        # Check if this oxygen is also attached to any atom in ring_j.
                        for nbr2 in nbr.GetNeighbors():
                            if nbr2.GetIdx() in ring_j:
                                connectivity += 1
                                found_connection = True
                                break
                    if found_connection:
                        break
                if found_connection:
                    break
    connected = (connectivity > 0)
    
    summary = (f"Detected {uronic_count} uronic acid-like and {glyco_count} glycosamine-like pyranose ring(s) "
               f"across {total_sugar_rings} candidate sugar ring(s).")
    if connected:
        summary += " Candidate rings appear connected via glycosidic bonds."
    else:
        summary += " No glycosidic connectivity among candidate rings detected."
    
    # Require at least four candidate sugar rings overall.
    if total_sugar_rings < 4:
        reason = f"Not enough sugar-like pyranose rings detected (found {total_sugar_rings}, need at least 4)."
        return False, reason
    
    # Require at least two rings of each type and roughly balanced counts.
    if uronic_count < 2 or glyco_count < 2:
        reason = (f"Repeating unit counts are too low: {uronic_count} uronic acid-like and "
                  f"{glyco_count} glycosamine-like rings (need at least 2 of each).")
        return False, reason
    
    if abs(uronic_count - glyco_count) > 1:
        reason = (f"Repeating unit counts are not balanced: {uronic_count} uronic acid-like vs. "
                  f"{glyco_count} glycosamine-like rings (difference > 1).")
        return False, reason
    
    # Mucopolysaccharides are typically large; here require a minimal molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 1000:
        return False, f"Molecular weight too low for a typical mucopolysaccharide ({mw:.1f} Da)."
    
    return True, summary

# Example usage (this test example is illustrative only):
if __name__ == "__main__":
    test_smiles = "O=C1NC(C(=O)NCCCC(=O)O)CC1"  # Replace with an actual mucopolysaccharide SMILES if available.
    result, reason = is_mucopolysaccharide(test_smiles)
    print("Result:", result)
    print("Reason:", reason)