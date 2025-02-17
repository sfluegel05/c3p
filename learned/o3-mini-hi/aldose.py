"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: Aldose – Aldehydic parent sugars (polyhydroxy aldehydes H[CH(OH)]nC(=O)H, n >= 2) and their intramolecular hemiacetals.
We now require that the overall molecule be of a size consistent with a monosaccharide (3 to 8 carbons).
For open-chain aldoses we require a terminal aldehyde group, at least 3 hydroxyl groups, and the overall carbon count in range.
For cyclic aldoses we require a 5- or 6-membered ring with exactly one ring oxygen, at least 2 exocyclic hydroxyl groups on ring carbons,
and again that the overall carbon count is small.
"""

from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose (either an open-chain polyhydroxy aldehyde or its cyclic hemiacetal form)
    based on its SMILES string.
    
    Improving on the previous version:
      - We first check that the molecule is not “too big” (typical parent sugars have 3-8 carbons).
      - For open-chain: we require the presence of a terminal aldehyde group ([CH1](=O))
        and at least three hydroxyl groups (as measured by [OX2H]).
      - For cyclic forms (typically furanose/pyranose rings): we require a 5- or 6-membered ring with exactly one oxygen,
        with the remaining atoms being carbons and with at least 2 exocyclic hydroxyl groups attached.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as an aldose, False otherwise.
        str: Reason for the classification.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total carbon atoms in the molecule
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(c_atoms)
    # Parent sugars (aldoses) are small molecules (typically 3-8 carbons)
    if not (3 <= n_carbons <= 8):
        return False, f"Total carbon count ({n_carbons}) is not in the typical range for a parent aldose (3-8)"
    
    # Count hydroxyl groups (-OH) using a SMARTS pattern [OX2H]
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    n_oh = len(oh_matches)
    
    # Try open-chain aldose: must have a terminal aldehyde and at least 3 hydroxyl groups.
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        if n_oh >= 3:
            return True, "Open-chain molecule with aldehyde and sufficient hydroxyl groups detected"
        else:
            return False, f"Open-chain aldehyde detected but only {n_oh} hydroxyl groups found (need at least 3)"
    
    # Otherwise, try cyclic (hemiacetal) structure.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Check each ring that is 5 or 6 atoms in size
    for ring in rings:
        if len(ring) not in (5, 6):
            continue
        # Count oxygens in the ring and collect the ring carbon indices
        ring_oxygen_count = 0
        ring_carbon_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                ring_carbon_indices.append(idx)
        
        # In a typical carbohydrate ring (pyranose or furanose) there is exactly one ring oxygen.
        if ring_oxygen_count != 1:
            continue
        
        # For a furanose ring, there are 4 carbons; for pyranose, 5 carbons.
        if len(ring_carbon_indices) not in (4, 5):
            continue
        
        # Check for exocyclic hydroxyl groups at the ring carbons.
        exo_oh_count = 0
        for idx in ring_carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            # Look at neighbors outside the ring
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue
                # Look for an -OH group connected: oxygen with at least one hydrogen
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    exo_oh_count += 1
                    break  # count only one per ring carbon
        
        if exo_oh_count >= 2:
            return True, "Cyclic hemiacetal structure consistent with an aldose detected (sugar ring with hydroxyls)"
    
    return False, "No open-chain aldehyde or typical cyclic monosaccharide (hemiacetal) pattern detected"

# For testing you might call:
# print(is_aldose("OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"))  # alpha-D-glucose (True)
# print(is_aldose("O=C1OC[C@@H]2O[C@@H](O[C@@H]3..."))  # should be rejected as too large