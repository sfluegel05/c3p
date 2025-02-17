"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: Aldose – Aldehydic parent sugars (polyhydroxy aldehydes H[CH(OH)]nC(=O)H, n >= 2) and their intramolecular hemiacetals.
"""

from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose (open-chain polyhydroxy aldehyde or its intramolecular hemiacetal)
    based on its SMILES string.
    
    The criteria are:
      - For open-chain aldoses: Must contain a terminal aldehyde group (–CHO), at least three hydroxyl (–OH) groups,
        and at least three carbon atoms.
      - For cyclic (hemiacetal) forms: Looks for a 5- or 6-membered ring that contains exactly one oxygen (the ring oxygen)
        and several hydroxyl groups attached to the ring carbons.
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule qualifies as an aldose, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total carbon atoms; an aldose must have at least 3 carbons.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 3:
        return False, "Not enough carbon atoms to be an aldose"
    
    # Count hydroxyl groups (-OH) using a simple SMARTS pattern.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(oh_pattern)
    if len(hydroxyl_matches) < 3:
        return False, "Not enough hydroxyl groups for a polyhydroxy aldehyde (need at least three -OH groups)"
    
    # First, check for an open-chain aldehyde.
    # The SMARTS pattern looks for a CHO group (a carbon with one hydrogen bound to a =O).
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Found an aldehyde. This meets the open-chain part of the definition.
        return True, "Open-chain molecule with aldehyde and sufficient hydroxyl groups detected"
    
    # If no free aldehyde is found, check for a cyclic structure (intramolecular hemiacetal).
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Look at each ring and try to find a typical carbohydrate ring
    for ring in atom_rings:
        ring_size = len(ring)
        # Only consider 5- or 6-membered rings (furanose or pyranose rings)
        if ring_size not in (5, 6):
            continue

        # Count oxygen atoms within the ring.
        ring_o_count = 0
        ring_carbon_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_o_count += 1
            elif atom.GetAtomicNum() == 6:
                ring_carbon_indices.append(idx)
        # In typical aldose rings there is exactly one ring oxygen.
        if ring_o_count != 1:
            continue

        # Now count hydroxyl groups attached to the ring carbons.
        ring_oh_count = 0
        for idx in ring_carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            # For each carbon, check neighbors for an -OH group.
            for neighbor in atom.GetNeighbors():
                # Exclude atoms that are also in the ring to consider only exocyclic hydroxyls.
                if neighbor.GetIdx() in ring:
                    continue
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    ring_oh_count += 1
                    break  # Count one -OH per ring carbon
        
        # It is common for a monosaccharide ring to have hydroxyl groups attached
        # to at least 2 (often more) of its ring carbons.
        if ring_oh_count >= 2:
            return True, "Cyclic hemiacetal structure consistent with an aldose detected (sugar ring with hydroxyls)"
    
    return False, "No open-chain aldehyde or typical cyclic monosaccharide (hemiacetal) pattern detected"