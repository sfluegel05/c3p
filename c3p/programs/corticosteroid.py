"""
Classifies: CHEBI:50858 corticosteroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_corticosteroid(smiles: str):
    """
    Determines if a molecule is a corticosteroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corticosteroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a steroid core structure (cyclopentanoperhydrophenanthrene nucleus)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Identify rings of size 6 and 5
    six_membered_rings = [ring for ring in atom_rings if len(ring) == 6]
    five_membered_rings = [ring for ring in atom_rings if len(ring) == 5]
    
    if len(six_membered_rings) < 3 or len(five_membered_rings) < 1:
        return False, "No steroid core structure found"

    # Check for the presence of functional groups typical of corticosteroids
    functional_groups = {
        "hydroxyl": Chem.MolFromSmarts("[OH]"),
        "carbonyl": Chem.MolFromSmarts("[C]=[O]"),
        "fluoro": Chem.MolFromSmarts("[F]"),
        "chloro": Chem.MolFromSmarts("[Cl]"),
        "acetate": Chem.MolFromSmarts("C(=O)O")
    }

    found_groups = {fg: mol.HasSubstructMatch(pattern) for fg, pattern in functional_groups.items()}

    if not found_groups["hydroxyl"] and not (found_groups["fluoro"] or found_groups["chloro"]):
        return False, "No hydroxyl or halogen group found"
    if not found_groups["carbonyl"]:
        return False, "No carbonyl group found"

    return True, "Molecule is classified as a corticosteroid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50858',
                          'name': 'corticosteroid',
                          'definition': 'A natural or synthetic analogue of '
                                        'the hormones secreted by the adrenal '
                                        'gland.',
                          'parents': ['CHEBI:35341']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 24,
    'num_false_positives': 10,
    'num_true_negatives': 10,
    'num_false_negatives': 3,
    'precision': 0.7058823529411765,
    'recall': 0.8888888888888888,
    'f1': 0.7868852459016393,
    'accuracy': None}