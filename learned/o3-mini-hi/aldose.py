"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: Aldose – Aldehydic parent sugars (polyhydroxy aldehydes H[CH(OH)]nC(=O)H, n>=2) and their 
intramolecular hemiacetals.

This implementation makes use of several rules:
1. The molecule must contain only H, C and O.
2. The number of carbon atoms must be between 4 and 8 and the overall count of oxygen atoms must equal
   the number of carbons (i.e. a CₙOₙ formula) – this removes deoxy-sugars and oxidized derivatives.
3. If the molecule is open‐chain (i.e. no suitable ring is found), then it must have a single terminal aldehyde
   (with the SMARTS “[CH1](=O)”) and (n – 1) –OH groups.
4. Otherwise, if at least one 5‐ or 6‐membered ring is found that contains exactly one oxygen and (importantly)
   none of its carbon atoms is in a carbonyl (C=O) environment, then we classify the molecule as a cyclic aldose.
5. If neither open‐chain nor cyclic criteria are met, the molecule is rejected as not a parent aldose.
   
Note: This is only one heuristic solution. Other implementations are possible.
"""
from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is a parent aldose based on its SMILES string.
    A valid parent aldose (open-chain or cyclic hemiacetal) according to our heuristic should:
      • only contain atoms H, C, and O;
      • have a total carbon count of 4 to 8;
      • have a molecular formula of CₙOₙ (i.e. oxygen count equals carbon count);
      • either (a) be an open-chain molecule with one terminal aldehyde group (SMARTS "[CH1](=O)") 
          and exactly (n-1) hydroxyl (-OH) groups, or 
          (b) contain a single 5- or 6-membered ring that has exactly one oxygen, and no ring carbon is 
          involved in a carbonyl bond.
          
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): Tuple of (True, reason) if classified as an aldose; (False, reason) if not.
                   If classification cannot be done, returns (None, None).
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Rule 1: Check that molecule contains only allowed elements (H, C, O)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in (1, 6, 8):
            return False, f"Molecule contains atom {atom.GetSymbol()}, not typical for a parent aldose"
    
    # Rule 2: Count total carbons and oxygens.
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Parent aldoses are typically 4- to 8-carbon sugars.
    if not (4 <= n_carbons <= 8):
        return False, f"Total carbon count ({n_carbons}) is not in the typical range for a parent aldose (4-8)"
    
    # Do a simple formula check: for a parent aldose, expect CₙOₙ.
    if n_oxygens != n_carbons:
        return False, f"Molecular formula has C{n_carbons}O{n_oxygens}, not typical for a parent aldose (should be CnOn)"
    
    # Prepare patterns
    aldehyde_pat = Chem.MolFromSmarts("[CH1](=O)")
    oh_pat = Chem.MolFromSmarts("[OX2H]")  # matches hydroxyl groups
    
    ring_info = mol.GetRingInfo()
    candidate_ring = None
    # Look for a 5- or 6-membered ring that might be a sugar ring.
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6):
            continue
        # Count oxygens in the ring and record carbon indices.
        ring_oxygen_count = 0
        ring_carbon_idx = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                ring_carbon_idx.append(idx)
        # A typical furanose or pyranose ring has exactly one oxygen.
        if ring_oxygen_count == 1:
            # Further require that none of the ring carbons is in a carbonyl (C=O) environment.
            carbonyl_in_ring = False
            for idx in ring_carbon_idx:
                atom = mol.GetAtomWithIdx(idx)
                for neighbor in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if neighbor.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                        carbonyl_in_ring = True
                        break
                if carbonyl_in_ring:
                    break
            if not carbonyl_in_ring:
                candidate_ring = ring
                break  # We found a qualifying cyclic sugar ring.
                
    # If we found a qualifying ring, classify as cyclic aldose.
    if candidate_ring is not None:
        return True, "Cyclic hemiacetal structure consistent with an aldose detected (sugar ring with hydroxyls)"
    
    # Otherwise, try to classify as an open-chain aldose.
    # For an open-chain parent aldose we expect:
    #   (a) a terminal aldehyde group; and 
    #   (b) every other carbon bears an -OH group.
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pat)
    if len(aldehyde_matches) != 1:
        # If no terminal aldehyde is found (or more than one), not an open-chain aldose.
        return False, f"Expected one terminal aldehyde group in open-chain aldose, found {len(aldehyde_matches)}"
    
    # Check that the aldehyde carbon is terminal.
    aldehyde_atom_idx = aldehyde_matches[0][0]  # the carbon that matches [CH1](=O)
    aldehyde_atom = mol.GetAtomWithIdx(aldehyde_atom_idx)
    # In an open-chain sugar the aldehyde carbon should have only one heavy-atom neighbor.
    if aldehyde_atom.GetDegree() != 1:
        return False, "Aldehyde group is not terminal"
    
    # Count hydroxyl groups.
    oh_matches = mol.GetSubstructMatches(oh_pat)
    # For an open-chain sugar with n carbons, we expect one aldehyde (no OH on that carbon)
    # and every other carbon (n - 1 total) must carry one OH.
    if len(oh_matches) != (n_carbons - 1):
        return False, f"Open-chain aldehyde detected but found {len(oh_matches)} hydroxyl groups (expected {n_carbons-1})"
    
    return True, "Open-chain molecule with terminal aldehyde and the expected number of hydroxyl groups detected"


# Example test calls:
if __name__ == '__main__':
    test_smiles = [
        # True positives (cyclic and open-chain aldoses):
        "O[C@H]1COC(O)[C@H](O)[C@@H]1O",            # L-lyxopyranose (cyclic)
        "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",     # D-allopyranose (cyclic)
        "OC[C@@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",  # alpha-L-gulose (cyclic)
        "OC[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@H]1O",   # alpha-D-gulose (cyclic)
        "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O",   # beta-D-gulose (cyclic)
        "[C@@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)CO",  # beta-L-idopyranose (cyclic)
        "[H]C([H])([C@]([H])(O)C=O)[C@]([H])(O)[C@@]([H])(C)O",  # tyvelose (open-chain)
        "OC[C@@H](O)[C@@H](O)C=O",                   # D-erythrose (open-chain)
        "O1[C@H](O)[C@H](O)[C@H](O)[C@@H](O)C1",      # beta-L-lyxopyranose (cyclic)
        "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",  # alpha-D-glucose (cyclic)
        "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO", # aldehydo-D-galactose (open-chain)
        # False positives:
        "OC(C(O)CC([H])=O)CO",                       # 2-deoxypentose (formula fails: C5 vs O4)
        "[H]C(=O)C(\\CC(O)=O)=C/C=C(\\O)C(O)=O",       # 5-formyl-2-hydroxyhepta-2,4-dienedioic acid (formula off)
        "[H][C@@]1(OC(=O)C(O)=C1O)[C@@H](O)CO",       # L-ascorbic acid (ring contains a carbonyl)
        # False negatives:
        "[H]C(=O)[C@H](O)CO",                        # D-glyceraldehyde (3-carbons, disqualified)
        "[H]C(=O)[C@@H](O)CO"                         # L-glyceraldehyde (3-carbons, disqualified)
    ]
    for s in test_smiles:
        result, reason = is_aldose(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*50}")