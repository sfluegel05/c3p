"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a fatty acid containing only methyl branches.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group at terminal position
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[OH]")
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if not matches:
        return False, "No carboxylic acid group found"

    # Exclude molecules with rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring structures"

    # Exclude molecules with functional groups other than carboxylic acid
    allowed_atoms = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Molecule contains heteroatoms not characteristic of fatty acids"

    # Check for functional groups other than carboxylic acid
    undesirable_groups = [
        Chem.MolFromSmarts("C(=O)N"),  # Amides
        Chem.MolFromSmarts("N"),       # Amines
        Chem.MolFromSmarts("O[H]"),    # Alcohols (exclude carboxylic OH)
        Chem.MolFromSmarts("[#6]=[O]"),# Ketones and aldehydes
    ]
    for group in undesirable_groups:
        if mol.HasSubstructMatch(group):
            return False, "Molecule contains functional groups other than carboxylic acid"

    # Identify the carboxyl carbon atom
    carboxyl_carbons = [match[0] for match in matches]  # First atom in pattern is the carbon
    carboxyl_c = carboxyl_carbons[0]

    # Find the longest carbon chain ending at the carboxyl carbon
    def find_longest_chain(atom_idx, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom_idx in visited or atom.GetAtomicNum() != 6:
            return []
        visited.add(atom_idx)
        paths = []
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 6:
                sub_path = find_longest_chain(n_idx, visited.copy())
                paths.append([atom_idx] + sub_path)
        if not paths:
            return [atom_idx]
        else:
            longest = max(paths, key=len)
            return longest

    longest_chain = find_longest_chain(carboxyl_c, set())
    main_chain_atoms = set(longest_chain)

    # Check for branches off the main chain
    for atom_idx in longest_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            continue
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in main_chain_atoms:
                branch = mol.GetAtomWithIdx(n_idx)
                # Check if branch is a methyl group
                if branch.GetAtomicNum() != 6:
                    return False, "Branch contains non-carbon atoms"
                # Methyl group should have only one heavy atom neighbor (the main chain carbon)
                heavy_neighbors = [nbr for nbr in branch.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if len(heavy_neighbors) != 1:
                    return False, "Branch larger than methyl group found"
                # Ensure no further branching from the methyl group
                for nbr in branch.GetNeighbors():
                    if nbr.GetIdx() != atom_idx and nbr.GetAtomicNum() != 1:
                        return False, "Branch larger than methyl group found"

    return True, "Molecule is a methyl-branched fatty acid with only methyl branches"

__metadata__ = {
    'chemical_class': {
        'name': 'methyl-branched fatty acid',
        'definition': 'Any branched-chain fatty acid containing methyl branches only.'
    },
    'message': None,
    'success': True,
    'error': '',
}