"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: aldose
"""

from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is an aldehydic parent sugar (polyhydroxy aldehyde) or its intramolecular hemiacetal form
    (a cyclic structure like a furanose or pyranose ring).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule contains only C, H, O atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 1, 8):
            return False, "Contains elements other than C, H, and O"

    # Check number of carbon atoms (aldoses typically have 3 to 7 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 2:
        return False, f"Too few carbons ({c_count}) for an aldose"

    # Define patterns
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")  # Aldehyde group
    hemiacetal_pattern = Chem.MolFromSmarts("[OX2H][CX4H1][OX2][CX4H][CX4H][OX2H]")  # Simplified hemiacetal pattern
    ring_5_or_6 = Chem.MolFromSmarts("[OX2H][C][C][C][O,C]")  # 5 or 6 membered ring with oxygen

    # Check for aldehyde group (open-chain form)
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Open-chain aldose
        # Check that all carbons except aldehyde carbon have hydroxyl groups
        aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
        aldehyde_c_idx = aldehyde_matches[0][0]
        aldehyde_c = mol.GetAtomWithIdx(aldehyde_c_idx)

        # Get chain carbons starting from aldehyde carbon
        visited = set()
        to_visit = [aldehyde_c]
        while to_visit:
            atom = to_visit.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() != 6:
                return False, "Chain contains non-carbon atoms"
            # For carbons other than aldehyde carbon, check for hydroxyl group
            if atom.GetIdx() != aldehyde_c_idx:
                hydroxyls = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() == 1]
                if len(hydroxyls) < 1:
                    return False, "A carbon in the chain does not have a hydroxyl group"

            # Add neighboring carbons
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                    to_visit.append(nbr)

        if len(visited) != c_count:
            return False, "Not all carbons are part of a single chain"
        return True, "Open-chain aldose with aldehyde group and hydroxylated carbons"

    else:
        # Check for cyclic hemiacetal forms (furanose or pyranose rings)
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 0:
            return False, "Does not contain aldehyde group or cyclic hemiacetal ring"

        # Find rings of size 5 or 6 containing one oxygen
        found_ring = False
        rings = ring_info.AtomRings()
        for ring in rings:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            ring_size = len(ring)
            if ring_size not in (5, 6):
                continue
            o_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
            if o_count != 1:
                continue
            # Check that the rest of the ring atoms are carbons
            if not all(atom.GetAtomicNum() == 6 for atom in ring_atoms if atom.GetAtomicNum() != 8):
                continue
            found_ring = True
            break  # Found suitable ring

        if not found_ring:
            return False, "No suitable furanose or pyranose ring found"

        # Check that most carbons have hydroxyl groups
        hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
        num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
        if num_hydroxyls < (c_count - 2):  # Allowing for anomeric carbon and possible deoxy sugars
            return False, f"Too few hydroxyl groups ({num_hydroxyls}) for an aldose"

        # Optionally check for exocyclic hydroxymethyl group
        # (common in aldoses forming pyranose rings)
        exocyclic_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
        valid_exocyclic = False
        for atom in exocyclic_carbons:
            # Check if carbon is CH2OH group
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 3:
                continue
            attached_to_ring = any(nbr.IsInRing() for nbr in neighbors if nbr.GetAtomicNum() == 6)
            hydroxyls = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() == 1]
            if attached_to_ring and len(hydroxyls) == 1:
                valid_exocyclic = True
                break

        # If no exocyclic CH2OH group, still acceptable
        return True, f"Cyclic aldose (hemiacetal form) with ring size {ring_size} and sufficient hydroxyl groups"

    return False, "Does not match aldose criteria"