"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: Prenols – any alcohol possessing the general formula 
    H-[CH2C(Me)=CHCH2]nOH 
where the carbon skeleton is exclusively composed of one or more isoprene units.
This version tries to “cover” the entire molecule (apart from hydrogens) with isoprene units.
Heuristic steps:
  1. Parse the SMILES string; the molecule must be acyclic.
  2. Only carbon and oxygen atoms are allowed.
  3. At least one free hydroxyl group (–OH not in –C(=O)OH) must be present.
     Moreover, at least one free –OH must be terminal (attached to a primary C).
  4. Count all carbon atoms; the total should be a multiple of five.
  5. Using two SMARTS patterns representing an isoprene unit:
         Pattern1: "[CH3]-[C]([CH3])=[CH]-[CH2]"
         Pattern2: "[CH2]-[C]([CH3])=[CH]-[CH2]"
     we gather all matches and then try to select the maximum number of non–overlapping matches.
  6. We “expect” that number to exactly equal (total_C / 5).
If everything checks, we classify the molecule as a prenol.
"""

from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    
    A prenol is defined as an alcohol whose carbon skeleton is exclusively
    composed of one or more isoprene units. (Isoprene unit is defined as 
    CH2–C(CH3)=CH–CH2.)
    
    Heuristic steps:
      - The molecule must be valid and acyclic.
      - It must contain only C and O as heavy atoms.
      - It must contain at least one free hydroxyl group (–OH) that is not part
        of a carboxyl.
      - It must have 1 or 2 free (i.e. terminal) hydroxyls.
      - Its total number of carbons should be divisible by 5.
      - It must “contain” n non–overlapping isoprene units (as determined by 
        two SMARTS patterns) where n = (number of carbons)//5.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule is classified as a prenol, False otherwise.
      str: Explanation of the decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject cyclic molecules.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic, not a prenol"
    
    # Reject if any heavy atom is not C (atomic num 6) or O (atomic num 8).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8):
            return False, f"Molecule contains atom {atom.GetSymbol()} not in (C,O)"
    
    # Identify all free hydroxyl groups.
    # We first match any hydroxyl: [OX2H]
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if not alcohol_matches:
        return False, "No hydroxyl (alcohol) group found"
    
    # Exclude those that are part of a carboxyl group: C(=O)[OX2H]
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX2H]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    # Convert matches to a set of oxygen atom indices that belong to carboxyl.
    carboxyl_oxygens = set()
    for match in carboxyl_matches:
        # match[1] is the oxygen in "O" of the pattern
        carboxyl_oxygens.add(match[1])
    
    # Filter free alcohols (oxygens not in a carboxyl)
    free_alcohols = [match for match in alcohol_matches if match[0] not in carboxyl_oxygens]
    if not free_alcohols:
        return False, "Only carboxylic acid hydroxyl(s) found; no free hydroxyl present"
    
    # Now ensure at least one of the free hydroxyls is “terminal”
    # i.e. the oxygen is attached to a primary carbon (the attached carbon has only one heavy neighbor aside from O).
    terminal_OH_count = 0
    for match in free_alcohols:
        # Each alcohol match gives the oxygen atom index.
        o_atom = mol.GetAtomWithIdx(match[0])
        # Get neighbors – there should be exactly one heavy neighbor.
        neighbors = [nbr for nbr in o_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]  # ignore hydrogens
        if len(neighbors) == 1:
            c_atom = neighbors[0]
            # Check how many heavy neighbors does the carbon have (excluding the oxygen we came from).
            c_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != o_atom.GetIdx()]
            # In a primary carbon (CH2 or CH3) attached to a chain, typically only one heavy neighbor.
            if len(c_neighbors) == 1:
                terminal_OH_count += 1
    if terminal_OH_count < 1:
        return False, "No terminal (primary) free hydroxyl group found"
    if terminal_OH_count > 2:
        return False, f"Too many terminal free hydroxyl groups found ({terminal_OH_count}); expected 1 or 2"
    
    # Count total number of carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_count = len(carbon_atoms)
    if carbon_count < 5:
        return False, f"Too few carbons ({carbon_count}) to be a prenol"
    if carbon_count % 5 != 0:
        return False, f"Total carbon count ({carbon_count}) is not a multiple of 5; does not match isoprene repeats"
    
    expected_units = carbon_count // 5  # number of isoprene units expected
    
    # Define two SMARTS for an isoprene unit.
    # Pattern1: CH3-C(CH3)=CH-CH2 (chain starts with CH3)
    isoprene_pat1 = Chem.MolFromSmarts("[CH3]-[C]([CH3])=[CH]-[CH2]")
    # Pattern2: CH2-C(CH3)=CH-CH2 (chain starts with CH2)
    isoprene_pat2 = Chem.MolFromSmarts("[CH2]-[C]([CH3])=[CH]-[CH2]")
    
    matches1 = mol.GetSubstructMatches(isoprene_pat1)
    matches2 = mol.GetSubstructMatches(isoprene_pat2)
    
    # Combine the matches (each match is a tuple of atom indices corresponding to the pattern)
    all_matches = list(matches1) + list(matches2)
    if not all_matches:
        return False, "No isoprene unit substructure match found"
    
    # To avoid counting overlapping isoprene units, we determine the maximum set of non-overlapping matches.
    # We use a simple recursive selection algorithm.
    def max_nonoverlapping(matches, used=set(), start=0):
        best = 0
        for i in range(start, len(matches)):
            # if this match overlaps with already used atoms, skip it.
            if set(matches[i]) & used:
                continue
            new_used = used | set(matches[i])
            count = 1 + max_nonoverlapping(matches, new_used, i+1)
            if count > best:
                best = count
        return best
    total_isoprene = max_nonoverlapping(all_matches)
    
    if total_isoprene != expected_units:
        return False, (f"Unexpected number of isoprene units: found {total_isoprene} non-overlapping match(es) "
                       f"but expected {expected_units} (from {carbon_count} carbons)")
    
    return True, (f"Classified as prenol: contains {expected_units} isoprene unit(s), "
                  f"{len(free_alcohols)} free hydroxyl group(s) (with {terminal_OH_count} terminal), "
                  f"and {carbon_count} carbons.")

# Example usage (for testing):
if __name__ == "__main__":
    # Use one of the examples (prenol itself: CC(C)=CCO)
    test_smiles = "CC(C)=CCO"
    result, reason = is_prenols(test_smiles)
    print("Result:", result)
    print("Reason:", reason)