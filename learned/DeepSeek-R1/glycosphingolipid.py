"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: glycosphingolipid (CHEBI:24404)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    Criteria: Carbohydrate attached via glycosidic linkage to ceramide/sphingoid base.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for presence of a carbohydrate (sugar) moiety
    # Look for a ring with at least 3 hydroxyl groups and an oxygen in the ring (pyranose/furanose)
    sugar_pattern = Chem.MolFromSmarts("[O;R]1@[C;R]@[C;R]@[C;R]@[C;R]@[C;R]1(-[OH])")  # pyranose
    if not mol.HasSubstructMatch(sugar_pattern):
        sugar_pattern_furan = Chem.MolFromSmarts("[O;R]1@[C;R]@[C;R]@[C;R]@[C;R]1(-[OH])")  # furanose
        if not mol.HasSubstructMatch(sugar_pattern_furan):
            return False, "No sugar moiety detected"

    # Check for ceramide/sphingoid base connected via O-glycosidic bond to sugar
    # Pattern: O (glycosidic) connected to sphingoid carbon and sugar's anomeric carbon
    # Sphingoid part: amino group (N) in a chain, ceramide has N-C=O
    glycosphingo_pattern = Chem.MolFromSmarts(
        "[NX3]-[CX4]-[CX4](-[OX2]-[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)"  # ceramide/sphingoid with O-sugar
    )
    if not mol.HasSubstructMatch(glycosphingo_pattern):
        # Alternative pattern for sphingoid (NH2 instead of amide)
        sphingoid_pattern = Chem.MolFromSmarts(
            "[NH2]-[CX4]-[CX4](-[OX2]-[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)"
        )
        if not mol.HasSubstructMatch(sphingoid_pattern):
            return False, "No glycosidic linkage between sugar and ceramide/sphingoid"

    # Verify long chain in ceramide/sphingoid (at least 12 carbons in fatty acid or sphingoid base)
    # Find amide groups (ceramide) or amine groups (sphingoid)
    amide_n = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and any(bond.GetBondType() == Chem.BondType.SINGLE and bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6 for bond in atom.GetBonds())]
    amine_n = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() >= 2]

    long_chain_found = False
    for n_atom in amide_n + amine_n:
        # Traverse the chain connected to N to check length
        # For ceramide: follow the amide's carbon chain
        # For sphingoid: follow the chain from N
        visited = set()
        queue = [(n_atom, 0)]
        max_chain_length = 0
        while queue:
            atom, depth = queue.pop(0)
            if atom in visited:
                continue
            visited.add(atom)
            if depth > max_chain_length:
                max_chain_length = depth
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor not in visited:
                    queue.append((neighbor, depth + 1))
        if max_chain_length >= 12:  # Arbitrary threshold for long chain
            long_chain_found = True
            break

    if not long_chain_found:
        return False, "Long chain in ceramide/sphingoid not found"

    return True, "Carbohydrate attached via glycosidic bond to ceramide/sphingoid"