"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    A carotenoid is a tetraterpenoid (C40), derived from psi,psi-carotene by various modifications
    such as hydrogenation, dehydrogenation, cyclization, oxidation, or combinations thereof.
    Retinoids are excluded.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude retinoids by checking for retinoid substructure
    retinoid_pattern = Chem.MolFromSmarts("CC(=O)C1=CC=CC=C1")
    if mol.HasSubstructMatch(retinoid_pattern):
        return False, "Molecule is a retinoid, which is excluded"

    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Molecule has {c_count} carbon atoms, expected at least 30 for carotenoids"

    # Detect long conjugated polyene chain (at least 10 alternating single and double bonds)
    pattern = Chem.MolFromSmarts("C(=C-C)*=C")
    matches = mol.GetSubstructMatches(pattern)
    if len(matches) < 10:
        return False, f"Found {len(matches)} conjugated double bonds, expected at least 10"

    # Check for beta-ionone ring or similar cyclic end groups
    beta_ionone_pattern = Chem.MolFromSmarts("C1(C=CC=CC=1)C")
    if not mol.HasSubstructMatch(beta_ionone_pattern):
        return False, "No beta-ionone ring found"

    # Allow for various functional groups (hydroxyl, keto, epoxide)
    # Do not exclude molecules lacking these groups

    # Check that molecule is derived from psi,psi-carotene skeleton
    # Using a simplified pattern for psi,psi-carotene
    psi_psi_carotene_pattern = Chem.MolFromSmarts("C(=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/CC)CC")
    if not mol.HasSubstructMatch(psi_psi_carotene_pattern):
        return False, "Molecule does not match psi,psi-carotene backbone"

    return True, "Molecule matches criteria for a carotenoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23044',
                          'name': 'carotenoid',
                          'definition': 'One of a class of tetraterpenoids (C40), formally derived from the acyclic parent, psi,psi-carotene by hydrogenation, dehydrogenation, cyclization, oxidation, or combination of these processes. This class includes carotenes, xanthophylls and certain compounds that arise from rearrangement of the skeleton of psi,psi-carotene or by loss of part of this structure. Retinoids are excluded.',
                          'parents': ['CHEBI:26859']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
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
    'accuracy': None }