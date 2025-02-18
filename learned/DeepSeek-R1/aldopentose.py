"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for explicit aldehyde group (O=CH-)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if aldehyde_matches:
        # Check each aldehyde group for being part of a 5-carbon chain
        for ald_atom_idx in [match[0] for match in aldehyde_matches]:
            # Get the longest chain starting from the aldehyde carbon
            # Assuming the aldehyde is at the end of the chain
            # Find the longest path from this atom excluding the aldehyde oxygen
            paths = Chem.GetAllPathsBetweenAtoms(mol, ald_atom_idx, maxPathLength=5)
            longest_chain = []
            for path in paths:
                if len(path) > len(longest_chain):
                    longest_chain = path
            # Check if the longest chain has 5 carbons (including aldehyde)
            if len(longest_chain) == 5:
                # Count hydroxyls on the chain (excluding aldehyde oxygen)
                hydroxyl_count = 0
                for atom_idx in longest_chain[1:]:  # exclude aldehyde carbon
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetAtomicNum() == 6:  # Carbon
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                                hydroxyl_count += 1
                                break
                if hydroxyl_count >= 3:
                    return True, "Contains aldehyde group and sufficient hydroxyls on 5-carbon chain"
        return False, "Aldehyde present but not part of a 5-carbon chain with sufficient hydroxyls"

    # Check for cyclic hemiacetal form (pentose sugar)
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        ring_size = len(ring)
        if ring_size not in (5, 6):
            continue
        oxygen_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxygen_in_ring != 1:
            continue

        # Check total carbons in molecule is 5
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if c_count != 5:
            continue

        # Count hydroxyl groups on ring carbons and hydroxymethyl group
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # Carbon in ring
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                        hydroxyl_count += 1
                        break  # Count once per carbon
        # Check for hydroxymethyl group (CH2OH) attached to the ring
        # This is the C5 in pentose
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6 and atom.GetDegree() == 2 and atom.GetTotalNumHs() == 2:
                # Possible hydroxymethyl group
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                        hydroxyl_count += 1
                        break

        if (ring_size == 5 and hydroxyl_count >= 3) or (ring_size == 6 and hydroxyl_count >= 4):
            return True, "Cyclic pentose with potential aldehyde group"

    return False, "Not an aldopentose"