"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: CHEBI: inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    Must have inositol linked via phosphodiester to a ceramide with amide and long chains.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Generalized inositol pattern (cyclohexane with multiple hydroxyls)
    inositol_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1(-[OH])O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol-like structure"

    # Find phosphate connected to inositol via phosphodiester
    # Pattern: inositol-O-P(=O)(O)-O-ceramide
    diester_pattern = Chem.MolFromSmarts("[C][O]P(=O)([O])[O][C]")
    if not mol.HasSubstructMatch(diester_pattern):
        return False, "No phosphodiester bridge"

    # Verify phosphate connects inositol to ceramide components
    # Get phosphate atoms
    phosphate_atoms = [match[1] for match in mol.GetSubstructMatches(Chem.MolFromSmarts("[O]P(=O)([O])[O]"))]
    bridge_valid = False
    for p_atom in phosphate_atoms:
        neighbors = [n for n in mol.GetAtomWithIdx(p_atom).GetNeighbors()]
        # Check two oxygen neighbors are connected to different parts (inositol and ceramide)
        inositol_connected = False
        ceramide_connected = False
        for n in neighbors:
            if n.GetAtomicNum() != 8:
                continue
            # Check if oxygen is connected to inositol
            for bond in n.GetBonds():
                other = bond.GetOtherAtom(n)
                if other.GetIdx() in [a.GetIdx() for a in mol.GetSubstructAtoms(inositol_pattern)]:
                    inositol_connected = True
            # Check if oxygen is connected to ceramide (amide with long chain)
            if not ceramide_connected:
                # Follow this oxygen to find amide group
                for b in n.GetBonds():
                    next_atom = b.GetOtherAtom(n)
                    if next_atom.GetAtomicNum() == 6:
                        amide_match = Chem.MolFromSmarts("[CX3](=O)[NX3H]")
                        if mol.GetSubstructMatch(amide_match, useQueryQuery=True):
                            # Check for long chains on both sides of amide
                            fatty_acid = next_atom
                            sphingoid_base = mol.GetAtomWithIdx(amide_match.GetMatches(mol)[0][2])
                            # Check chain lengths (at least 12 carbons each)
                            def chain_length(atom):
                                visited = set()
                                stack = [(atom, 0)]
                                max_length = 0
                                while stack:
                                    a, depth = stack.pop()
                                    if a.GetIdx() in visited or a.GetAtomicNum() != 6:
                                        continue
                                    visited.add(a.GetIdx())
                                    max_length = max(max_length, depth)
                                    for neighbor in a.GetNeighbors():
                                        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != a.GetIdx():
                                            stack.append((neighbor, depth + 1))
                                return max_length
                            if chain_length(fatty_acid) >= 12 and chain_length(sphingoid_base) >= 12:
                                ceramide_connected = True
        if inositol_connected and ceramide_connected:
            bridge_valid = True
            break
    if not bridge_valid:
        return False, "Phosphodiester doesn't link inositol to ceramide"

    # Check for hydroxyl in sphingoid base (adjacent to amide nitrogen)
    amide_n = [match[2] for match in mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[NX3H]"))]
    hydroxyl_found = False
    for n_idx in amide_n:
        n_atom = mol.GetAtomWithIdx(n_idx)
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() > 0:
                hydroxyl_found = True
                break
        if not hydroxyl_found:
            # Check within 3 bonds of nitrogen
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, 3, n_idx)
            atoms = set([n_idx])
            for bond_idx in env:
                atoms.add(mol.GetBondWithIdx(bond_idx).GetBeginAtomIdx())
                atoms.add(mol.GetBondWithIdx(bond_idx).GetEndAtomIdx())
            for a_idx in atoms:
                if mol.GetAtomWithIdx(a_idx).GetAtomicNum() == 8 and mol.GetAtomWithIdx(a_idx).GetTotalNumHs() > 0:
                    hydroxyl_found = True
                    break
        if hydroxyl_found:
            break
    if not hydroxyl_found:
        return False, "No hydroxyl in sphingoid base"

    return True, "Inositol-phosphoceramide with phosphodiester bridge and ceramide features"