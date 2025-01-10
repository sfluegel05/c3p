"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: CHEBI:28053 aldopentose
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a monosaccharide with five carbon atoms and an aldehyde group at one end,
    which can exist in open-chain or cyclic forms (furanose or pyranose rings).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldopentose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove explicit hydrogens for consistency
    mol = Chem.RemoveHs(mol)

    # Count number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 5:
        return False, f"Contains {num_carbons} carbon atoms, should be 5 for an aldopentose"

    # Exclude molecules with disallowed functional groups
    disallowed_patterns = [
        '[CX3](=O)[OX1H0-,OX2H1]',  # Carboxylic acids and derivatives
        'P(=O)(O)(O)',              # Phosphate groups
        'C1OC(=O)[C@@H](O)[C@H]1O',  # Lactone rings
    ]
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Contains disallowed functional groups (acid, phosphate, lactone)"

    # Check for open-chain aldopentose with terminal aldehyde group
    open_chain_pattern = Chem.MolFromSmarts('[O=CH][CH2][CH](O)[CH](O)CO')
    if mol.HasSubstructMatch(open_chain_pattern):
        return True, "Contains open-chain aldopentose structure with terminal aldehyde group"

    # Check for furanose ring (5-membered ring with oxygen)
    furanose_pattern = Chem.MolFromSmarts('C1OC([CH](O))[CH](O)[CH2]O1')
    if mol.HasSubstructMatch(furanose_pattern):
        return True, "Contains furanose ring aldopentose"

    # Check for pyranose ring (6-membered ring with oxygen)
    pyranose_pattern = Chem.MolFromSmarts('C1OC([CH](O))[CH](O)[CH](O)[CH2]O1')
    if mol.HasSubstructMatch(pyranose_pattern):
        return True, "Contains pyranose ring aldopentose"

    # Check for hemiacetal linkage in cyclic forms
    hemiacetal_pattern = Chem.MolFromSmarts('C1([OH])OC([CH](O))[CH](O)[CH](O)O1')
    if mol.HasSubstructMatch(hemiacetal_pattern):
        return True, "Contains cyclic aldopentose with hemiacetal linkage"

    # If none of the patterns matched, it's not an aldopentose
    return False, "Does not match aldopentose structures"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:28053',
        'name': 'aldopentose',
        'definition': 'A pentose with a (potential) aldehyde group at one end.',
        'parents': ['CHEBI:16842', 'CHEBI:24869']
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