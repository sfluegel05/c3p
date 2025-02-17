"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: Aldose – Aldehydic parent sugars (polyhydroxy aldehydes H[CH(OH)]nC(=O)H, n>=2) and their 
intramolecular hemiacetals.

Improvements over the previous version:
  1. Ensures that the molecule is a small parent sugar (3 to 8 carbons).
  2. Rejects molecules with atoms other than H, C, and O (to avoid cases with phosphate, sulfate, selenium, etc).
  3. For open-chain aldoses, requires a terminal aldehyde group and at least three hydroxyl groups.
  4. For cyclic aldoses, requires a 5- or 6-membered ring containing exactly one ring oxygen (the rest carbons)
     plus at least two exocyclic hydroxyl groups attached to ring carbons.
"""

from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose (either an open‐chain polyhydroxy aldehyde or its cyclic
    hemiacetal form) based on its SMILES string.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule qualifies as an aldose, False otherwise.
      str: Reason for the classification result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Rule 1: Check that molecule contains only allowed elements (H, C, O).
    for atom in mol.GetAtoms():
        Z = atom.GetAtomicNum()
        if Z not in (1, 6, 8):  # only allow H, C, O
            return False, f"Molecule contains atom {atom.GetSymbol()}, not typical for a parent aldose"
    
    # Rule 2: Count total carbon atoms, which should be from 3 to 8 for a parent sugar.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(carbon_atoms)
    if not (3 <= n_carbons <= 8):
        return False, f"Total carbon count ({n_carbons}) is not in the typical range for a parent aldose (3-8)"
    
    # Count hydroxyl groups: we use the SMARTS for an -OH group.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    n_oh = len(oh_matches)
    
    # Rule 3: Try open-chain aldose: must have a terminal aldehyde and at least three hydroxyl groups.
    # We use a SMARTS pattern for an aldehyde where the carbon bears one hydrogen.
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        if n_oh >= 3:
            return True, "Open-chain molecule with terminal aldehyde and sufficient hydroxyl groups detected"
        else:
            return False, f"Open-chain aldehyde detected but only {n_oh} hydroxyl groups found (need at least 3)"
    
    # Rule 4: Otherwise, try cyclic (hemiacetal) aldose.
    # Look through all rings in the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    for ring in rings:
        # Consider only 5- or 6-membered rings
        if len(ring) not in (5, 6):
            continue
        
        # Count how many oxygens are present in the ring.
        ring_oxygen_count = 0
        ring_carbon_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                ring_carbon_indices.append(idx)
        # For a typical furanose/pyranose ring, expect exactly one ring oxygen.
        if ring_oxygen_count != 1:
            continue
        
        # Optionally, we could check that the number of ring carbons is (ring size - 1)
        if len(ring_carbon_indices) != (len(ring) - 1):
            continue
        
        # Count exocyclic hydroxyl groups on ring carbons.
        exo_oh_count = 0
        for idx in ring_carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                # Only consider neighbors outside the ring.
                if neighbor.GetIdx() in ring:
                    continue
                # Look for an -OH group (oxygen with at least one hydrogen attached)
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    exo_oh_count += 1
                    break  # count maximum one per ring carbon
        if exo_oh_count >= 2:
            return True, "Cyclic hemiacetal structure consistent with an aldose detected (sugar ring with hydroxyls)"
    
    return False, "No open-chain aldehyde or typical cyclic monosaccharide (hemiacetal) pattern detected"


# Example test calls:
if __name__ == '__main__':
    test_smiles = [
        # True positives (cyclic/alcohol sugars):
        "O[C@H]1COC(O)[C@H](O)[C@@H]1O",            # L-lyxopyranose
        "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",     # D-allopyranose
        "OC[C@@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",  # alpha-L-gulose
        # Open-chain example:
        "OC[C@@H](O)[C@@H](O)C=O",                   # D-erythrose
        # False positive candidate (should be rejected because of extra patterns)
        "OC(C(O)CC([H])=O)CO",                       # 2-deoxypentose example (likely too deoxygenated)
    ]
    for s in test_smiles:
        result, reason = is_aldose(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*40}")