"""
Classifies: CHEBI:23066 cephalosporin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for a 4-membered beta-lactam ring
    beta_lactam_ring = None
    for ring in rings.AtomRings():
        if len(ring) == 4:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetSymbol() == 'C' or atom.GetSymbol() == 'N' for atom in atoms):
                beta_lactam_ring = ring
                break

    if beta_lactam_ring is None:
        return False, "No 4-membered beta-lactam ring found"

    # Check for a 6-membered side ring attached to the beta-lactam ring
    side_ring = None
    for ring in rings.AtomRings():
        if len(ring) == 6 and any(atom_idx in beta_lactam_ring for atom_idx in ring):
            side_ring = ring
            break

    if side_ring is None:
        return False, "No 6-membered side ring attached to the beta-lactam ring found"

    # Check if the 6-membered side ring contains sulfur
    if not any(mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'S' for atom_idx in side_ring):
        return False, "6-membered side ring does not contain sulfur"

    return True, "Molecule is a cephalosporin"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23066',
                          'name': 'cephalosporin',
                          'definition': 'A class of beta-lactam antibiotics '
                                        'differing from the penicillins in '
                                        'having a 6-membered, rather than a '
                                        '5-membered, side ring.  Although '
                                        'cephalosporins are among the most '
                                        'commonly used antibiotics in the '
                                        'treatment of routine infections, and '
                                        'their use is increasing over time, '
                                        'they can cause a range of '
                                        'hypersensitivity reactions, from '
                                        'mild, delayed-onset cutaneous '
                                        'reactions to life-threatening '
                                        'anaphylaxis in patients with '
                                        'immunoglobulin E (IgE)-mediated '
                                        'allergy.',
                          'parents': ['CHEBI:38311']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 11,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}