"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a pentose (5 carbons in the main chain) with a (potential) aldehyde group.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for explicit aldehyde group (O=CH-)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Verify presence of at least 3 hydroxyl groups (characteristic of aldopentoses)
        hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
        if hydroxyl_count >= 3:
            return True, "Contains aldehyde group and sufficient hydroxyls"
        else:
            return False, "Aldehyde present but insufficient hydroxyl groups"

    # Check for cyclic hemiacetal form (pentose sugar)
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        ring_size = len(ring)
        # Check for pentose-typical ring sizes (furanose:5, pyranose:6) with exactly one oxygen
        if ring_size not in (5, 6):
            continue
        oxygen_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxygen_in_ring != 1:
            continue  # Require exactly one oxygen (hemiacetal)

        # Count hydroxyl groups on ring carbons
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # Carbon in ring
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                        hydroxyl_count += 1
                        break  # Count once per carbon

        # Validate hydroxyl count based on ring size
        if (ring_size == 5 and hydroxyl_count >= 3) or (ring_size == 6 and hydroxyl_count >= 4):
            return True, "Cyclic pentose with potential aldehyde group"

    return False, "Not an aldopentose"