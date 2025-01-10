"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:25099 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    An icosanoid is any signaling molecule arising from oxidation of the three C20 essential fatty acids
    (EFAs): icosapentaenoic acid (EPA), arachidonic acid (AA), and dihomo-gamma-linolenic acid (DGLA).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 22:
        return False, f"Molecule has {c_count} carbon atoms; expected approximately 20 for an icosanoid"

    # Count the number of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Molecule has insufficient oxygen atoms; expected oxidation products"

    # Check for presence of functional groups
    functional_groups = []

    # Hydroxyl group (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        functional_groups.append("hydroxyl")

    # Ketone group (>C=O)
    ketone_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if mol.HasSubstructMatch(ketone_pattern):
        functional_groups.append("ketone")

    # Epoxide group (three-membered ring with an oxygen)
    epoxide_pattern = Chem.MolFromSmarts("[C]1[O][C]1")
    if mol.HasSubstructMatch(epoxide_pattern):
        functional_groups.append("epoxide")

    # Peroxide group (-O-O-)
    peroxide_pattern = Chem.MolFromSmarts("[OX2][OX2]")
    if mol.HasSubstructMatch(peroxide_pattern):
        functional_groups.append("peroxide")

    if not functional_groups:
        return False, "Molecule lacks characteristic functional groups of icosanoids (hydroxyl, ketone, epoxide, peroxide)"

    # Estimate the number of rotatable bonds to infer chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, f"Not enough rotatable bonds ({n_rotatable}); molecule may not have long aliphatic chains"

    # Check for cyclopentane ring (for prostaglandins)
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")
    has_cyclopentane = mol.HasSubstructMatch(cyclopentane_pattern)

    # Check for conjugated double bonds (for leukotrienes)
    conjugated_diene_pattern = Chem.MolFromSmarts("C=CC=CC=CC=C")
    has_conjugated_diene = mol.HasSubstructMatch(conjugated_diene_pattern)

    if not (has_cyclopentane or has_conjugated_diene):
        return False, "Molecule lacks characteristic rings or conjugated double bonds of icosanoids"

    return True, f"Molecule has characteristics of an icosanoid (carbon count: {c_count}, functional groups: {', '.join(functional_groups)})"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25099',
        'name': 'icosanoid',
        'definition': 'Any member of the group of signalling molecules arising from oxidation of the three C20 essential fatty acids (EFAs) icosapentaenoic acid (EPA), arachidonic acid (AA) and dihomo-gamma-linolenic acid (DGLA).',
        'parents': []
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
    'attempt': 1,
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