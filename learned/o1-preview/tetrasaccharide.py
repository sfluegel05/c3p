"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for pyranose and furanose rings (generic monosaccharide units)
    # These patterns match 5 or 6-membered rings with oxygen and hydroxyl groups
    pyranose_pattern = Chem.MolFromSmarts('C1[C@@H]([O])[C@H](O)[C@@H](O)[C@H](O)O1')
    furanose_pattern = Chem.MolFromSmarts('C1[C@H](O)[C@@H](O)[C@H](O)O1')

    # Find monosaccharide units
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    total_monosaccharides = len(pyranose_matches) + len(furanose_matches)

    if total_monosaccharides != 4:
        return False, f"Found {total_monosaccharides} monosaccharide units, expected 4"

    # Collect atoms involved in monosaccharide units
    mono_atom_indices = set()
    for match in pyranose_matches + furanose_matches:
        mono_atom_indices.update(match)

    # Identify glycosidic linkages
    glycosidic_bonds = []
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        idx1 = atom1.GetIdx()
        idx2 = atom2.GetIdx()

        # Check if bond is between two monosaccharide units via an oxygen atom
        if idx1 in mono_atom_indices and idx2 in mono_atom_indices:
            # Exclude bonds within the same monosaccharide unit
            in_same_unit = False
            for match in pyranose_matches + furanose_matches:
                if idx1 in match and idx2 in match:
                    in_same_unit = True
                    break
            if in_same_unit:
                continue

            # Check if bond is an oxygen bridge
            if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or \
               (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8):
                glycosidic_bonds.append(bond)

    if len(glycosidic_bonds) != 3:
        return False, f"Found {len(glycosidic_bonds)} glycosidic linkages, expected 3"

    # Check connectivity between monosaccharide units via glycosidic linkages
    from collections import defaultdict, deque

    # Build a mapping from atom index to monosaccharide unit
    atom_to_unit = {}
    for unit_idx, match in enumerate(pyranose_matches + furanose_matches):
        for idx in match:
            atom_to_unit[idx] = unit_idx

    # Build connectivity graph
    graph = defaultdict(set)
    for bond in glycosidic_bonds:
        idx1 = bond.GetBeginAtom().GetIdx()
        idx2 = bond.GetEndAtom().GetIdx()
        unit1 = atom_to_unit.get(idx1)
        unit2 = atom_to_unit.get(idx2)
        if unit1 is not None and unit2 is not None and unit1 != unit2:
            graph[unit1].add(unit2)
            graph[unit2].add(unit1)

    # Check if the monosaccharide units are connected
    visited = set()
    queue = deque()
    queue.append(0)
    visited.add(0)
    while queue:
        current = queue.popleft()
        for neighbor in graph[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    if len(visited) != 4:
        return False, "Monosaccharide units are not properly connected via glycosidic linkages"

    return True, "Molecule contains 4 monosaccharide units connected via glycosidic linkages"