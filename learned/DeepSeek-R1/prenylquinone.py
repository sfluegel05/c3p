"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: CHEBI: prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    Requires a 1,4-benzoquinone or 1,4-naphthoquinone core with at least one
    polyprenyl-derived side chain containing multiple isoprene units.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Corrected quinone core patterns allowing substitutions
    # 1,4-benzoquinone: two ketones in para positions on benzene ring
    benzoquinone = Chem.MolFromSmarts("[#6]1(-[#8]=[#6]):[#6](-[#8]=[#6]):[#6]:[#6]:[#6]:[#6]:1")
    # 1,4-naphthoquinone: two ketones in para positions on naphthalene
    naphthoquinone = Chem.MolFromSmarts("[#6]12:[#6](-[#8]=[#6]):[#6]:[#6]:[#6]:[#6]:[#6]:1:[#6](-[#8]=[#6]):[#6]:[#6]:2")

    # Check for quinone core existence
    has_quinone = mol.HasSubstructMatch(benzoquinone) or mol.HasSubstructMatch(naphthoquinone)
    if not has_quinone:
        return False, "No quinone core detected"

    # Flexible prenyl chain pattern: at least two isoprene units with branching
    # Matches (CH3)-C-CH2-C-CH2- pattern with possible double bonds
    prenyl_pattern = Chem.MolFromSmarts(
        "[CH3]-[CH2](-[CH3])-[CH2]-[CH2](-[CH3])"
    )

    # Find all possible attachment points on quinone core
    core_matches = mol.GetSubstructMatches(benzoquinone) or mol.GetSubstructMatches(naphthoquinone)

    # Check all possible prenyl chain attachments (direct or via oxygen)
    for core_match in core_matches:
        core_atoms = set(core_match)
        for atom_idx in core_atoms:
            core_atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in core_atom.GetNeighbors():
                # Follow oxygen links (e.g., methoxy attachments)
                if neighbor.GetAtomicNum() == 8:
                    for oxy_neighbor in neighbor.GetNeighbors():
                        if oxy_neighbor.GetAtomicNum() == 6:
                            chain_root = oxy_neighbor
                            break
                    else:
                        continue
                elif neighbor.GetAtomicNum() == 6:
                    chain_root = neighbor
                else:
                    continue

                # Traverse chain and check for prenyl characteristics
                visited = set()
                stack = [(chain_root, 0)]
                max_chain_length = 0
                methyl_count = 0

                while stack:
                    atom, depth = stack.pop()
                    if atom.GetIdx() in visited:
                        continue
                    visited.add(atom.GetIdx())

                    # Track chain length and methyl branches
                    if depth > max_chain_length:
                        max_chain_length = depth
                    if any(n.GetAtomicNum() == 6 and n.GetTotalNumHs() >= 3 for n in atom.GetNeighbors()):
                        methyl_count += 1

                    # Continue traversal
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() not in core_atoms and nbr.GetAtomicNum() == 6:
                            stack.append((nbr, depth + 1))

                # Require minimum chain characteristics
                if max_chain_length >= 8 and methyl_count >= 2:
                    return True, "Quinone core with polyprenyl side chain"

    # Check for conjugated isoprene systems using double bond count
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds >= 3 and rdMolDescriptors.CalcNumAmideBonds(mol) == 0:
        return True, "Multiple conjugated systems typical of prenyl chains"

    return False, "No polyprenyl chain attached to quinone core"