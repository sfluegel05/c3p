"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: Sterol Ester
Definition: A steroid ester obtained by formal condensation of the carboxy group 
of any carboxylic acid with the 3-hydroxy group of a sterol.

Heuristic:
  1. The molecule must contain an ester functional group (a C(=O)O fragment).
  2. The molecule must have a steroid nucleus—that is, at least one 5‐membered ring 
     and at least three 6‐membered rings.
  3. Among the ester groups, identify the “ester oxygen” (the alcohol side, i.e. the oxygen attached
     to the carbonyl carbon by a single bond). Check that one of its neighbors (other than the carbonyl carbon)
     is in one of the rings from the steroid nucleus.
     
This improved heuristic should help remedy cases in which the ester oxygen is not directly in a ring.
"""
from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    
    A sterol ester is defined as a steroid ester formed by reaction of any 
    carboxylic acid with the 3-hydroxy group of a sterol. In this case, we expect
    an ester (C(=O)O) present and that the oxygen from the ester (which used to be the hydroxyl from the sterol)
    is attached (directly or via a short bond) to a steroid nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a sterol ester, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, enforce the presence of an ester group.
    # We use a SMARTS pattern for an ester: a carbon atom with a double bond to oxygen,
    # and a single bond to an oxygen atom.
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional group (C(=O)O) found"
    
    # Next, verify the presence of a steroid nucleus.
    # A typical sterol has a cyclopentanoperhydrophenanthrene nucleus (four fused rings),
    # which entails at least one 5-membered ring and three 6-membered rings.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    num_5_membered = sum(1 for ring in rings if len(ring) == 5)
    num_6_membered = sum(1 for ring in rings if len(ring) == 6)
    if num_5_membered < 1 or num_6_membered < 3:
        return False, ("Steroid nucleus not detected – expected at least one 5-membered ring and three 6-membered rings, "
                       f"found {num_5_membered} and {num_6_membered}, respectively")
    
    # Build a set of atoms in 5- or 6-membered rings (heuristically, steroid rings).
    steroid_atoms = set()
    for ring in rings:
        if len(ring) in (5, 6):
            steroid_atoms.update(ring)
    
    # Now, among the ester groups we search for one that likely comes from the condensation
    # at the 3-hydroxy of a sterol.
    # To do this, for each ester match, we:
    #   - Identify the carbonyl carbon (atomic number 6) in the match.
    #   - Among the oxygen atoms in the match, pick the one that is attached via a single bond 
    #     to the carbonyl carbon (this is the ester oxygen from the alcohol part).
    #   - Then check that one of its neighbors (other than the carbonyl carbon) is in the steroid nucleus.
    sterol_ester_found = False
    for match in ester_matches:
        # match will contain the indices of atoms that match the SMARTS pattern.
        # Identify the candidate carbonyl carbon and the oxygen(s) in this match.
        candidate_c = None
        candidate_oxygens = []
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                candidate_c = idx
            elif atom.GetAtomicNum() == 8:
                candidate_oxygens.append(idx)
        # Skip if we did not find a clear carbon
        if candidate_c is None or not candidate_oxygens:
            continue
        
        # For each oxygen in candidate_oxygens, check the bond order to the carbon.
        # We want the oxygen that is connected by a single bond.
        for oxy_idx in candidate_oxygens:
            bond = mol.GetBondBetweenAtoms(candidate_c, oxy_idx)
            if bond is None:
                continue
            if bond.GetBondTypeAsDouble() != 1.0:  # only accept a single bond (not the carbonyl oxygen)
                continue
            ester_oxygen_idx = oxy_idx
            ester_oxygen = mol.GetAtomWithIdx(ester_oxygen_idx)
            # Get neighbors of the ester oxygen excluding the carbonyl carbon.
            neighbors = [nbr.GetIdx() for nbr in ester_oxygen.GetNeighbors() if nbr.GetIdx() != candidate_c]
            if not neighbors:
                continue
            # Check if any neighbor is part of the steroid nucleus
            # (i.e. the atom is in steroid_atoms)
            if any(nbr in steroid_atoms for nbr in neighbors):
                sterol_ester_found = True
                break
        if sterol_ester_found:
            break
    
    if not sterol_ester_found:
        return False, ("Ester group found but the oxygen from the ester does not appear to connect (directly) to any "
                       "atom in a steroid ring system (expected for a sterol-derived alcohol)")
    
    return True, "Contains an ester group with an alcohol side linked to a steroid ring system, consistent with a sterol ester"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided example SMILES (e.g., CE(20:2(6Z,9Z)))
    smiles_example = "O([C@@H]1CC=2[C@@](C3C(C4[C@@](C(CC4)[C@@H](CCCC(C)C)C)(CC3)C)CC2)(CC1)C)C(CCCCCCCCC/C=C\\C/C=C\\CCCCC)=O"
    result, reason = is_sterol_ester(smiles_example)
    print("Result:", result)
    print("Reason:", reason)