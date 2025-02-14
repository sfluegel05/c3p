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
    An aldose is an aldehydic sugar (polyhydroxy aldehyde) or its intramolecular hemiacetal form
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

    # Check that the molecule contains only C, H, O atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 1, 8):
            return False, "Contains elements other than C, H, and O"

    # Check the number of carbon atoms (aldoses typically have 3 to 7 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3 or c_count > 7:
        return False, f"Number of carbons ({c_count}) not typical for an aldose"

    # Check for other rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:
        return False, "Contains more than one ring"

    # Count the number of hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    num_hydroxyl = len(mol.GetSubstructMatches(hydroxyl_pattern))

    # Check for aldehyde group (open-chain form)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

    if aldehyde_matches:
        # Open-chain form
        # Each non-aldehyde carbon should have exactly one hydroxyl group
        # Get atoms connected to aldehyde carbon
        aldehyde_c_idx = aldehyde_matches[0][0]
        aldehyde_c = mol.GetAtomWithIdx(aldehyde_c_idx)
        neighbor_carbons = [nbr.GetIdx() for nbr in aldehyde_c.GetNeighbors() if nbr.GetAtomicNum() == 6]

        chain_carbons = neighbor_carbons.copy()
        visited = set([aldehyde_c_idx])
        while neighbor_carbons:
            current_idx = neighbor_carbons.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            if atom.GetAtomicNum() != 6:
                return False, "Chain contains non-carbon atoms"
            # Check that carbon has exactly one hydroxyl group
            hydroxyls = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() == 1]
            if len(hydroxyls) != 1:
                return False, "A carbon in the chain does not have exactly one hydroxyl group"
            # Add neighboring carbons
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                    neighbor_carbons.append(nbr.GetIdx())
        if len(visited) != c_count:
            return False, "Not all carbons are part of a single chain"
        return True, f"Open-chain aldose with aldehyde group and {num_hydroxyl} hydroxyl groups"
    else:
        # Check for cyclic hemiacetal forms (furanose or pyranose rings)
        if ring_info.NumRings() == 0:
            return False, "Does not contain aldehyde group or cyclic hemiacetal ring"

        # Only one ring allowed
        if ring_info.NumRings() > 1:
            return False, "Contains more than one ring"

        # Get ring atoms
        ring_atoms = ring_info.AtomRings()[0]
        ring_size = len(ring_atoms)
        if ring_size not in (5, 6):
            return False, f"Ring size {ring_size} is not typical for furanose or pyranose"

        # Check that the ring contains exactly one oxygen
        ring_oxygen_atoms = [idx for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(ring_oxygen_atoms) != 1:
            return False, f"Ring contains {len(ring_oxygen_atoms)} oxygen atoms instead of 1"

        # Check that all ring carbons have exactly one hydroxyl group
        for idx in ring_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                hydroxyls = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() == 1]
                if len(hydroxyls) != 1:
                    return False, "A ring carbon does not have exactly one hydroxyl group"

        # Check for exocyclic -CH2OH group (common in aldoses)
        exocyclic_atoms = set(mol.GetAtoms()) - set(mol.GetAtomWithIdx(idx) for idx in ring_atoms)
        exocyclic_carbons = [atom for atom in exocyclic_atoms if atom.GetAtomicNum() == 6]
        if len(exocyclic_carbons) > 1:
            return False, "Too many exocyclic carbons"
        elif len(exocyclic_carbons) == 1:
            # Check that exocyclic carbon is a -CH2OH group
            exo_c = exocyclic_carbons[0]
            neighbors = exo_c.GetNeighbors()
            if len(neighbors) != 3:
                return False, "Exocyclic carbon does not have 3 neighbors"
            hydroxyls = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() == 1]
            if len(hydroxyls) != 1:
                return False, "Exocyclic carbon does not have a hydroxyl group"
            attached_to_ring = any(nbr.GetIdx() in ring_atoms for nbr in neighbors if nbr.GetAtomicNum() == 6)
            if not attached_to_ring:
                return False, "Exocyclic carbon is not attached to the ring"
        elif len(exocyclic_carbons) == 0:
            pass  # No exocyclic carbons, acceptable
        else:
            return False, "Unexpected number of exocyclic carbons"

        return True, f"Cyclic aldose (hemiacetal form) with ring size {ring_size} and {num_hydroxyl} hydroxyl groups"

    return False, "Does not match aldose criteria"