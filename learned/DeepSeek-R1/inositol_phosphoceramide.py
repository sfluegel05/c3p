"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: CHEBI: inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    An inositol phosphoceramide consists of an inositol residue linked via a phosphodiester
    bridge to a ceramide moiety, which contains a sphingoid base and a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for inositol (myo-inositol core with multiple hydroxyls)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol moiety detected"

    # Check for phosphate group connected to inositol
    phosphate_pattern = Chem.MolFromSmarts("[O][P](=O)([O])[O]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Verify phosphate bridges inositol and ceramide
    # Find phosphate atoms and check connectivity
    p_connected_to_inositol = False
    p_connected_to_ceramide = False
    inositol_atoms = set(mol.GetSubstructMatch(inositol_pattern))
    for p_match in phosphate_matches:
        p_idx = p_match[0]
        p_atom = mol.GetAtomWithIdx(p_idx)
        neighbors = p_atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                # Check if oxygen is connected to inositol
                for bond in neighbor.GetBonds():
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetIdx() in inositol_atoms:
                        p_connected_to_inositol = True
                # Check if oxygen leads to ceramide (amide group)
                # Traverse from this oxygen to find amide
                visited = set()
                stack = [(neighbor, 0)]
                while stack:
                    current_atom, depth = stack.pop()
                    if current_atom.GetIdx() in visited:
                        continue
                    visited.add(current_atom.GetIdx())
                    if depth > 6:  # Limit search depth
                        continue
                    # Look for amide group (CONH)
                    for bond in current_atom.GetBonds():
                        other_atom = bond.GetOtherAtom(current_atom)
                        if other_atom.GetAtomicNum() == 6:  # Carbon
                            for b in other_atom.GetBonds():
                                if b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(other_atom).GetAtomicNum() == 8:
                                    # Check adjacent to NH
                                    for nb in other_atom.GetNeighbors():
                                        if nb.GetAtomicNum() == 7:  # Nitrogen
                                            p_connected_to_ceramide = True
                                            break
                    if not p_connected_to_ceramide:
                        stack.extend([(a, depth+1) for a in current_atom.GetNeighbors() if a.GetAtomicNum() != 15])  # Exclude P
        if p_connected_to_inositol and p_connected_to_ceramide:
            break
    else:
        return False, "Phosphate not bridging inositol and ceramide"

    # Check for ceramide structure (amide with long chains)
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3H]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group (ceramide) found"

    # Check fatty acid chain length (at least 14 carbons)
    fatty_acid_valid = False
    for amide in amide_matches:
        carbon_idx = amide[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        # Traverse away from amide to find fatty acid chain
        chain = []
        stack = [(carbon_atom, 0)]
        visited = set()
        while stack:
            atom, depth = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() != 6:
                continue
            chain.append(atom)
            # Stop if chain branches or hits heteroatom
            if depth > 20:  # Prevent infinite loops
                break
            for bond in atom.GetBonds():
                next_atom = bond.GetOtherAtom(atom)
                if next_atom.GetAtomicNum() == 6 and next_atom.GetIdx() not in visited:
                    stack.append((next_atom, depth + 1))
        if len(chain) >= 14:
            fatty_acid_valid = True
            break
    if not fatty_acid_valid:
        return False, "Fatty acid chain too short"

    # Check sphingoid base has hydroxyl group
    hydroxyl_found = False
    for amide in amide_matches:
        nitrogen_idx = amide[2]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        # Traverse from nitrogen to find hydroxyl
        stack = [(nitrogen_atom, 0)]
        visited = set()
        while stack:
            atom, depth = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if depth > 10:
                break
            if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1:  # Hydroxyl oxygen
                hydroxyl_found = True
                break
            for bond in atom.GetBonds():
                next_atom = bond.GetOtherAtom(atom)
                if next_atom.GetAtomicNum() in [6,8] and next_atom.GetIdx() not in visited:
                    stack.append((next_atom, depth + 1))
        if hydroxyl_found:
            break
    if not hydroxyl_found:
        return False, "Sphingoid base lacks hydroxyl group"

    return True, "Contains inositol-phosphoceramide structure with phosphodiester bridge, amide, and long chains"