"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: pyrimidine deoxyribonucleoside
"""

from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside consists of a pyrimidine base attached to a deoxyribose sugar via a β-N1-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove salts and small fragments
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, default=mol, key=lambda m: m.GetNumAtoms())

    # Exclude molecules with phosphate groups (nucleotides)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group, likely a nucleotide"

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    pyrimidine_ring_atoms = None
    # Identify pyrimidine rings (aromatic, six-membered, with exactly two nitrogen atoms)
    for ring in atom_rings:
        if len(ring) == 6:
            aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
            num_nitrogens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if aromatic and num_nitrogens == 2:
                pyrimidine_ring_atoms = set(ring)
                break

    if pyrimidine_ring_atoms is None:
        return False, "No pyrimidine base found"

    sugar_ring_atoms = None
    # Identify sugar rings (five-membered rings containing an oxygen atom)
    for ring in atom_rings:
        if len(ring) == 5:
            has_oxygen = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)
            if has_oxygen:
                sugar_ring_atoms = set(ring)
                break

    if sugar_ring_atoms is None:
        return False, "No sugar ring found"

    # Check for glycosidic bond between pyrimidine base nitrogen and sugar carbon
    glycosidic_bond_found = False
    for atom_idx in pyrimidine_ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 7:  # Nitrogen atom in base
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in sugar_ring_atoms and neighbor.GetAtomicNum() == 6:
                    # Check if this carbon is the anomeric carbon (attached to oxygen)
                    for nb in neighbor.GetNeighbors():
                        if nb.GetIdx() in sugar_ring_atoms and nb.GetAtomicNum() == 8:
                            glycosidic_bond_found = True
                            break
                if glycosidic_bond_found:
                    break
        if glycosidic_bond_found:
            break

    if not glycosidic_bond_found:
        return False, "No glycosidic bond between base nitrogen and sugar carbon"

    # Check that the sugar is deoxyribose (lacks 2' hydroxyl group)
    # In deoxyribose, the 2' carbon (next to ring oxygen) should not have an -OH group
    ring_oxygen_idx = None
    for idx in sugar_ring_atoms:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
            ring_oxygen_idx = idx
            break

    if ring_oxygen_idx is None:
        return False, "Sugar ring does not contain ring oxygen"

    # Find 2' carbon (attached to ring oxygen)
    two_prime_carbons = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(ring_oxygen_idx).GetNeighbors()
                         if nbr.GetIdx() in sugar_ring_atoms and nbr.GetAtomicNum() == 6]

    if len(two_prime_carbons) != 2:
        return False, "Sugar ring structure invalid"

    # Identify the 2' carbon (the one not connected to the glycosidic bond)
    for c_idx in two_prime_carbons:
        if c_idx not in [neighbor_idx,]:  # Exclude anomeric carbon
            two_prime_carbon_idx = c_idx
            break

    # Check if 2' carbon has an -OH group (should not in deoxyribose)
    two_prime_carbon = mol.GetAtomWithIdx(two_prime_carbon_idx)
    has_2prime_OH = False
    for nbr in two_prime_carbon.GetNeighbors():
        if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in sugar_ring_atoms:
            has_2prime_OH = True
            break

    if has_2prime_OH:
        return False, "Sugar is ribose (has 2' OH group), not deoxyribose"

    return True, "Contains pyrimidine base attached to deoxyribose via N-glycosidic bond"