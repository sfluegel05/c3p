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
    Must have inositol linked via phosphodiester to a ceramide with amide and long chains.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Generalized inositol (cyclohexanehexol or similar)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol"

    # Find phosphate connected to inositol
    phosphate = Chem.MolFromSmarts("[O]P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate)
    if not phosphate_matches:
        return False, "No phosphate"

    # Check phosphate bridges inositol and ceramide
    inositol_atoms = set(mol.GetSubstructMatch(inositol_pattern))
    bridge_found = False
    for p_match in phosphate_matches:
        p_idx = p_match[0]
        p_atom = mol.GetAtomWithIdx(p_idx)
        # Check connections: one O to inositol, another O to ceramide
        linked_to_inositol = False
        linked_to_ceramide = False
        for neighbor in p_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 8:
                continue
            # Check if neighbor is connected to inositol
            for bond in neighbor.GetBonds():
                other = bond.GetOtherAtom(neighbor)
                if other.GetIdx() in inositol_atoms:
                    linked_to_inositol = True
            # Check if neighbor connects to ceramide (amide with long chain)
            stack = [(neighbor, 0)]
            visited = set()
            while stack:
                current, depth = stack.pop()
                if current.GetIdx() in visited:
                    continue
                visited.add(current.GetIdx())
                if depth > 10:
                    break
                # Look for amide group (CONH)
                if current.GetAtomicNum() == 6 and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in current.GetBonds()):
                    for nbr in current.GetNeighbors():
                        if nbr.GetAtomicNum() == 7 and any(bond.GetBondType() == Chem.BondType.SINGLE for bond in nbr.GetBonds()):
                            # Check for long chain
                            chain_length = 0
                            chain_stack = [(current, 0)]
                            chain_visited = set()
                            while chain_stack:
                                catom, clevel = chain_stack.pop()
                                if catom.GetIdx() in chain_visited:
                                    continue
                                chain_visited.add(catom.GetIdx())
                                if catom.GetAtomicNum() == 6:
                                    chain_length += 1
                                    if chain_length >= 14:
                                        linked_to_ceramide = True
                                        break
                                    for b in catom.GetBonds():
                                        next_atom = b.GetOtherAtom(catom)
                                        if next_atom.GetAtomicNum() == 6 and clevel < 30:
                                            chain_stack.append((next_atom, clevel+1))
                                if linked_to_ceramide:
                                    break
                if not linked_to_ceramide:
                    for b in current.GetBonds():
                        next_a = b.GetOtherAtom(current)
                        if next_a.GetAtomicNum() in [6,8] and next_a.GetIdx() not in visited:
                            stack.append((next_a, depth+1))
        if linked_to_inositol and linked_to_ceramide:
            bridge_found = True
            break
    if not bridge_found:
        return False, "Phosphate bridge missing"

    # Verify ceramide has sphingoid base with hydroxyl
    amide_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[NX3H]"))
    if not amide_matches:
        return False, "No amide"
    hydroxyl_found = False
    for amide in amide_matches:
        nh_atom = mol.GetAtomWithIdx(amide[2])
        # Traverse from N to find OH
        stack = [(nh_atom, 0)]
        visited = set()
        while stack:
            atom, depth = stack.pop()
            if atom.GetIdx() in visited or depth > 6:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0:
                hydroxyl_found = True
                break
            for bond in atom.GetBonds():
                next_a = bond.GetOtherAtom(atom)
                if next_a.GetAtomicNum() in [6,8]:
                    stack.append((next_a, depth+1))
        if hydroxyl_found:
            break
    if not hydroxyl_found:
        return False, "No hydroxyl in sphingoid"

    return True, "Inositol-phosphoceramide with phosphodiester bridge and ceramide"