"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: Mucopolysaccharide (glycosaminoglycan)
Definition: “Any of the group of polysaccharides composed of alternating units from uronic acids
and glycosamines, and commonly partially esterified with sulfuric acid.”
Heuristic improvements:
  1. Parse molecule.
  2. Iterate over rings (only 5- or 6-membered rings) and require that each candidate sugar ring be “clean”
     in that all its atoms are either carbon (6) or oxygen (8) with exactly one ring oxygen.
  3. For each such ring, examine substituents on its ring carbons (neighbors not in the ring):
       a. If a substituent is nitrogen, mark the ring as glycosamine–like.
       b. If a substituent is a carbon that bears at least one double‐bonded oxygen plus at least one additional oxygen 
          (mimicking a –C(=O)O– group), then mark the ring as uronic acid–like.
     If a given ring shows both features (or is ambiguous), we skip that ring.
  4. Also search for sulfate ester groups with SMARTS.
  5. Finally, require at least two uronic and two glycosamine–like rings with nearly balanced counts (difference ≤ 1).
Note: This heuristic is still approximate.
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide (glycosaminoglycan) by seeking sugar-like rings 
    with substituents diagnostic of either uronic acid or glycosamine units. Also counts sulfate ester groups.
    
    Args:
      smiles (str): SMILES string of the molecule.

    Returns:
      bool: True if the molecule is classified as a mucopolysaccharide, False otherwise.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Counters for sugar units
    uronic_count = 0
    glyco_count = 0

    # Get all ring sets (each ring as a tuple of atom indices)
    rings = mol.GetRingInfo().AtomRings()

    # Loop through each ring
    for ring in rings:
        # Only consider typical sugar rings (5 or 6 atoms)
        if len(ring) not in (5, 6):
            continue
        
        # Get the atoms in this ring.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Check that every atom in the ring is either carbon (6) or oxygen (8)
        if any(atom.GetAtomicNum() not in (6, 8) for atom in ring_atoms):
            continue
        
        # A canonical sugar ring (pyranose or furanose) contains exactly one ring oxygen.
        ring_oxy = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        if ring_oxy != 1:
            continue
        
        # Initialize flags for this ring.
        uronic_flag = False
        glyco_flag = False
        
        # For each atom in the ring that is carbon, examine substituents not in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # only process ring carbons
            
            for nbr in atom.GetNeighbors():
                # Skip atoms that are in the ring
                if nbr.GetIdx() in ring:
                    continue
                # Check for glycosamine indicator: directly attached nitrogen
                if nbr.GetAtomicNum() == 7:
                    glyco_flag = True
                
                # Check for uronic acid indicator: substituent carbon with carboxyl group pattern.
                if nbr.GetAtomicNum() == 6:
                    # For the substituent carbon, look at its neighbors (excluding the ring atom)
                    ext_neighbors = [n for n in nbr.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
                    doubleO = 0
                    countO = 0
                    for ext in ext_neighbors:
                        if ext.GetAtomicNum() == 8:
                            countO += 1
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), ext.GetIdx())
                            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                                doubleO += 1
                    # Classical carboxyl group has at least one double-bonded oxygen and (total oxygens) at least 2.
                    if doubleO >= 1 and countO >= 2:
                        uronic_flag = True
        
        # If the ring shows both features, it is ambiguous; skip it.
        if uronic_flag and glyco_flag:
            continue
        
        # Increase counts as appropriate.
        if uronic_flag:
            uronic_count += 1
        elif glyco_flag:
            glyco_count += 1

    # Also, look for sulfate ester groups (a common modification)
    sulfate_smarts = "OS(=O)(=O)O"
    sulfate_query = Chem.MolFromSmarts(sulfate_smarts)
    sulfate_matches = mol.GetSubstructMatches(sulfate_query) if sulfate_query is not None else []
    n_sulfate = len(sulfate_matches)
    
    summary = (f"Detected {uronic_count} uronic acid-like ring(s) and "
               f"{glyco_count} glycosamine-like ring(s).")
    if n_sulfate > 0:
        summary += f" Also found {n_sulfate} sulfate ester group(s)."
    
    # Check that we have at least 2 of each and that the counts are roughly balanced (difference ≤ 1).
    if uronic_count < 2 or glyco_count < 2:
        reason = (f"Not enough sugar-like units: found {uronic_count} uronic acid-like and "
                  f"{glyco_count} glycosamine-like rings (need at least 2 of each)")
        return False, reason
    
    if abs(uronic_count - glyco_count) > 1:
        reason = (f"Repeating unit counts are not balanced (difference > 1): "
                  f"{uronic_count} uronic acid-like vs {glyco_count} glycosamine-like rings")
        return False, reason

    return True, summary

# Example usage:
if __name__ == "__main__":
    # Test with an example SMILES (this one is just illustrative).
    test_smiles = "O=C1NC(C(=O)NCCCC(=O)O)CC1"  # not a real mucopolysaccharide; replace as needed
    result, reason = is_mucopolysaccharide(test_smiles)
    print("Result:", result)
    print("Reason:", reason)