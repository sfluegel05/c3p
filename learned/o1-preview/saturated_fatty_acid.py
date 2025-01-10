"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: CHEBI:15841 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid is a fatty acid containing no carbon-to-carbon multiple bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH or deprotonated form -C(=O)O-)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O,H,-1]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Expected exactly one carboxylic acid group, found {len(carboxylic_acid_matches)}"

    # Check for carbon-carbon double bonds
    if mol.HasSubstructMatch(Chem.MolFromSmarts('C=C')):
        return False, "Contains carbon-carbon double bonds (unsaturation)"

    # Check for carbon-carbon triple bonds
    if mol.HasSubstructMatch(Chem.MolFromSmarts('C#C')):
        return False, "Contains carbon-carbon triple bonds (unsaturation)"

    # Check for ring structures
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains ring structures, not a fatty acid"

    # Check for heteroatoms other than O in carboxylic acid group
    allowed_atoms = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in allowed_atoms:
            return False, f"Contains heteroatoms other than C, H, and O (found atomic number {atomic_num})"
        # Exclude oxygen atoms not part of the carboxylic acid group
        if atomic_num == 8:
            atom_idx = atom.GetIdx()
            is_in_carboxyl = any(atom_idx in match for match in carboxylic_acid_matches)
            if not is_in_carboxyl:
                return False, "Contains oxygen atoms outside carboxylic acid group"

    # Check for additional functional groups beyond carboxylic acid
    # Exclude hydroxyl groups beyond carboxylic acid hydroxyl
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    # Count hydroxyl groups not in carboxylic acid
    num_hydroxyls = 0
    for match in hydroxyl_matches:
        atom_idx = match[0]
        is_in_carboxyl = any(atom_idx in match for match in carboxylic_acid_matches)
        if not is_in_carboxyl:
            num_hydroxyls += 1
    if num_hydroxyls > 0:
        return False, f"Contains {num_hydroxyls} additional hydroxyl group(s)"

    # Exclude amine groups
    amine_pattern = Chem.MolFromSmarts('[NX3;H2,H1,H0;!$(NC=O)]')
    if mol.HasSubstructMatch(amine_pattern):
        return False, "Contains amine group"

    # Exclude other functional groups like ketones, aldehydes, esters, amides, halogens
    forbidden_patterns = [
        Chem.MolFromSmarts('C=O'),        # Carbonyl (ketones, aldehydes)
        Chem.MolFromSmarts('C(=O)O[C,H]'),  # Esters
        Chem.MolFromSmarts('C(=O)N'),     # Amides
        Chem.MolFromSmarts('[F,Cl,Br,I]'), # Halogens
        Chem.MolFromSmarts('O=CN'),       # Amides
        Chem.MolFromSmarts('S'),          # Sulfur-containing groups
        Chem.MolFromSmarts('P'),          # Phosphorus-containing groups
    ]
    for pattern in forbidden_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains forbidden functional groups"

    # Get the carbon chain length (number of carbons excluding the carboxyl carbon)
    carboxyl_c_idx = carboxylic_acid_matches[0][0]
    visited = set()
    chain_length = get_aliphatic_chain_length(mol, carboxyl_c_idx, visited)
    if chain_length < 1:
        return False, f"Carbon chain too short for fatty acid (found {chain_length} carbon)"

    return True, "Molecule is a saturated fatty acid"

def get_aliphatic_chain_length(mol, atom_idx, visited):
    """
    Recursively counts the number of carbon atoms connected in the aliphatic chain.

    Args:
        mol: RDKit molecule object
        atom_idx: Index of the current atom
        visited: Set of visited atom indices

    Returns:
        int: Length of the carbon chain excluding the carboxyl carbon
    """
    visited.add(atom_idx)
    atom = mol.GetAtomWithIdx(atom_idx)
    length = 0
    for neighbor in atom.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx in visited:
            continue
        neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
        if neighbor_atom.GetAtomicNum() == 6:  # Carbon atom
            bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
            # Only consider single bonds
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            length += 1 + get_aliphatic_chain_length(mol, neighbor_idx, visited)
    return length

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:15841',
        'name': 'saturated fatty acid',
        'definition': 'Any fatty acid containing no carbon to carbon multiple bonds.',
        'parents': ['CHEBI:35366', 'CHEBI:26559']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}