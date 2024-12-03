"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside (Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid structure (four fused rings A/B/C/D)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not (ring_sizes.count(6) >= 3 and ring_sizes.count(5) >= 1):
        return False, "Molecule does not contain the required steroid structure"

    # Check for lactone ring
    lactone_found = False
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in atoms):
                lactone_found = True
                break

    if not lactone_found:
        return False, "Molecule does not contain a lactone ring"

    # Check for sugar residues (pyranose or furanose rings)
    sugar_found = False
    for ring in ring_info.AtomRings():
        if len(ring) in [5, 6]:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in atoms):
                if all(atom.GetSymbol() in ['C', 'O'] for atom in atoms):
                    sugar_found = True
                    break

    if not sugar_found:
        return False, "Molecule does not contain sugar residues"

    return True, "Molecule is classified as a cardiac glycoside"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83970',
                          'name': 'cardiac glycoside',
                          'definition': 'Steroid lactones containing sugar '
                                        'residues that act on the contractile '
                                        'force of the cardiac muscles.',
                          'parents': ['CHEBI:24400', 'CHEBI:26766']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}