"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: Aldose – Aldehydic parent sugars (polyhydroxy aldehydes H[CH(OH)]nC(=O)H, n>=2)
and their intramolecular hemiacetals.

This version contains improved heuristics:
 - Only allow H, C and O; require 4–8 carbons.
 - Reject any molecule with a ketone-type C=O (i.e. a C(=O) where the carbon has no hydrogen).
 - For an open–chain candidate, look for a terminal aldehyde defined by the SMARTS
   "[CX3H1](=O)" and confirm that (aside from the bonded oxygen) that carbon has exactly one heavy neighbor.
   Then count free hydroxyls ([OX2H]) and allow a tolerance of ±1 compared to (n_carbons – 1).
 - For a cyclic candidate, search for a 5– or 6–membered ring that has exactly one ring oxygen and no carbonyl on any ring carbon.
   Then count how many ring carbons have an exocyclic –OH (as detected by bonds to an O with H) and accept if that count
   is within one of the fully hydroxylated value.
 - If neither pattern fits, the molecule is rejected.
Note: Sugar chemistry is complex; this heuristic–based method will not catch every nuance.
"""

from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose sugar (either open-chain
    (with terminal aldehyde) or cyclic hemiacetal) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if classified as an aldose, 
                     False otherwise; the second element gives a reason.
                     If the classification cannot be attempted, (None, None) is returned.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Rule 1: Only allow atoms H, C and O.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, f"Molecule contains {atom.GetSymbol()}, which is not allowed for a parent aldose"
    
    # Count carbons – expect a parent aldose to have 4-8.
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (4 <= n_carbons <= 8):
        return False, f"Carbon count ({n_carbons}) is not in the expected 4-8 range for a parent aldose"
    
    # Define SMARTS for terminal aldehyde: the carbon must have exactly one hydrogen.
    aldehyde_smarts = "[CX3H1](=O)"
    aldehyde_query = Chem.MolFromSmarts(aldehyde_smarts)
    
    # SMARTS for free hydroxyl group.
    oh_smarts = "[OX2H]"
    oh_query = Chem.MolFromSmarts(oh_smarts)
    
    # Rule 2: Reject if any C=O (carbonyl) is detected on a carbon without a hydrogen
    # (i.e. likely part of a ketone or non-aldehydic carbonyl).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond.GetBondTypeAsDouble() == 2.0:
                        # if the carbonyl carbon does not have exactly one H, then not an aldehyde
                        if atom.GetTotalNumHs() != 1:
                            return False, "Non-terminal carbonyl (ketone) group detected, not an aldose"
    
    # Count free –OH groups in the molecule.
    oh_matches = mol.GetSubstructMatches(oh_query)
    n_oh = len(oh_matches)
    
    # For fully hydroxylated aldose (no deoxy substitution) we expect (n_carbons – 1) – but we allow ±1 tolerance.
    full_oh = n_carbons - 1

    # FIRST: Try the open-chain candidate.
    open_chain_matches = mol.GetSubstructMatches(aldehyde_query)
    open_chain_valid = False
    if len(open_chain_matches) >= 1:
        # We look for a terminal aldehyde: the carbon (first index of the match) should have exactly one heavy neighbor
        # apart from the carbonyl oxygen.
        for match in open_chain_matches:
            aldehyde_idx = match[0]
            aldehyde_atom = mol.GetAtomWithIdx(aldehyde_idx)
            # Count neighbors excluding oxygen that is double-bonded (the carbonyl oxygen)
            heavy_neighbors = [nbr for nbr in aldehyde_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
            # Remove the oxygen that is part of the C=O bond:
            heavy_neighbors = [nbr for nbr in heavy_neighbors 
                               if not (nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(aldehyde_idx, nbr.GetIdx()).GetBondTypeAsDouble() == 2.0)]
            if len(heavy_neighbors) == 1:
                open_chain_valid = True
                break
    # If we find a valid open-chain aldehyde and the –OH count is close to full hydroxylation, accept.
    if open_chain_valid:
        if abs(n_oh - full_oh) <= 1:
            return True, ("Open-chain aldehyde detected (terminal carbonyl) with free OH count " +
                          f"of {n_oh} (expected ~{full_oh})")
        else:
            reason = f"Open-chain aldehyde detected but free OH count {n_oh} is not within tolerance of expected {full_oh}"
            # we do not immediately reject—it might be a cyclic candidate instead.
    # SECOND: Try the cyclic (hemiacetal) candidate.
    ring_info = mol.GetRingInfo()
    candidate_ring = None
    candidate_ring_type = None
    for ring in ring_info.AtomRings():
        # Consider only 5- or 6-membered rings.
        if len(ring) not in (5, 6):
            continue
        # Count how many atoms in the ring are oxygen.
        ring_oxygen_count = 0
        ring_carbon_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                ring_carbon_indices.append(idx)
        # In a typical pyranose or furanose ring for an aldose, there is exactly one ring oxygen.
        if ring_oxygen_count != 1:
            continue
        
        # For each of the ring carbons, ensure none is involved in a carbonyl double bond.
        ring_ok = True
        for idx in ring_carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                    ring_ok = False
                    break
            if not ring_ok:
                break
        if not ring_ok:
            continue
        
        # Now, count how many ring carbons have an exocyclic -OH attached.
        oh_count_on_ring = 0
        for idx in ring_carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            # Loop over neighbors that are not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check if the neighbor oxygen is part of an -OH group.
                    if nbr.HasSubstructMatch(oh_query):
                        oh_count_on_ring += 1
                        break
        # For a full aldose ring, one expects each ring carbon to have an exocyclic OH.
        if abs(oh_count_on_ring - len(ring_carbon_indices)) <= 1:
            candidate_ring = ring
            candidate_ring_type = "cyclic"
            break

    if candidate_ring is not None:
        return True, ("Cyclic hemiacetal structure detected in a ring of size " +
                      f"{len(candidate_ring)} with {oh_count_on_ring} exocyclic OH groups (expected ~{len(candidate_ring)-1} on carbons)")
    
    # If neither open-chain nor cyclic criteria were fully met then reject.
    if open_chain_valid:
        return False, f"Open-chain candidate detected but free OH count {n_oh} is not within tolerance of expected {full_oh}"
    return False, "No valid open-chain aldehyde or cyclic hemiacetal pattern detected for a parent aldose"


# Example test calls:
if __name__ == '__main__':
    test_smiles = [
        # True positives (aldoses)
        "O[C@H]1COC(O)[C@H](O)[C@@H]1O",             # L-lyxopyranose (cyclic)
        "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",      # D-allopyranose (cyclic)
        "OC[C@@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",# alpha-L-gulose (cyclic)
        "OC[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@H]1O",    # alpha-D-gulose (cyclic)
        "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O",    # beta-D-gulose (cyclic)
        "[C@@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)CO",# beta-L-idopyranose (cyclic)
        "O1[C@H](O)[C@H](O)[C@H](O)[C@@H](O)C1",         # beta-L-lyxopyranose (cyclic)
        "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",     # alpha-D-glucose (cyclic)
        # Open chain candidates:
        "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",   # aldehydo-D-galactose (open-chain)
        "OC[C@@H](O)[C@@H](O)C=O",                      # D-erythrose (open-chain)
        # The following examples were previously false positives (ketoses, etc.) and should be rejected:
        "OC[C@H]1O[C@](O)(CO)[C@@H](O)[C@H](O)[C@@H]1O", # beta-D-sedoheptulopyranose (ketose)
        "OC[C@H]1OC(O)(CO)[C@@H](O)[C@H]1O",             # D-tagatofuranose (ketose)
    ]
    for s in test_smiles:
        result, reason = is_aldose(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*50}")