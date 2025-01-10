"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: CHEBI:52214 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is any fatty acid containing at least one aldehydic or ketonic group
    in addition to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Identify atoms in carboxylic acid groups to exclude
    carboxylic_acid_atoms = set()
    for match in carboxylic_acid_matches:
        carboxylic_acid_atoms.update(match)

    # Check for aldehyde groups (excluding carboxylic acid)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = []
    for match in mol.GetSubstructMatches(aldehyde_pattern):
        # Ensure the aldehyde is not part of the carboxylic acid
        if not any(atom_idx in carboxylic_acid_atoms for atom_idx in match):
            aldehyde_matches.append(match)

    # Check for ketone groups (excluding carboxylic acid)
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6,#1]")
    ketone_matches = []
    for match in mol.GetSubstructMatches(ketone_pattern):
        atom_indices = set(match)
        # Ensure the ketone is not part of the carboxylic acid or other excluded groups
        if not atom_indices & carboxylic_acid_atoms:
            carbonyl_c = mol.GetAtomWithIdx(match[1])
            # Exclude esters, amides by checking neighboring oxygens or nitrogens
            neighbors = carbonyl_c.GetNeighbors()
            is_simple_ketone = True
            for neigh in neighbors:
                if neigh.GetIdx() != match[0] and neigh.GetIdx() != match[2]:
                    if neigh.GetAtomicNum() in [7,8]:
                        is_simple_ketone = False
                        break
            if is_simple_ketone:
                ketone_matches.append(match)

    # If no aldehyde or ketone groups are found, it's not an oxo fatty acid
    if not aldehyde_matches and not ketone_matches:
        return False, "No aldehyde or ketone group found in addition to carboxylic acid"

    # Count the number of carbon atoms to verify it's a fatty acid
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 4:
        return False, f"Too few carbon atoms ({num_carbons}) for a fatty acid"

    # Check for other functional groups not typical in fatty acids
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            return False, "Contains nitrogen, not typical for fatty acids"
        if atom.GetAtomicNum() == 15:
            return False, "Contains phosphorus, not typical for fatty acids"
        if atom.GetAtomicNum() == 16:
            return False, "Contains sulfur, not typical for fatty acids"

    # Check for number of carboxylic acid groups (should be one)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Contains {len(carboxylic_acid_matches)} carboxylic acid groups, expected 1"

    # Calculate fraction of sp3 carbons (should be high in fatty acids)
    num_sp3_carbons = sum(1 for atom in mol.GetAtoms()
                          if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)
    sp3_fraction = num_sp3_carbons / num_carbons
    if sp3_fraction < 0.5:
        return False, f"Low fraction of sp3 carbons ({sp3_fraction:.2f}), molecule may not be aliphatic enough for fatty acid"

    # If all checks pass, it's an oxo fatty acid
    return True, "Contains carboxylic acid and aldehyde/ketone groups with sufficient aliphatic carbon chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:52214',
        'name': 'oxo fatty acid',
        'definition': 'Any fatty acid containing at least one aldehydic or ketonic group in addition to the carboxylic acid group.',
        'parents': ['CHEBI:25360', 'CHEBI:25352']
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
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}