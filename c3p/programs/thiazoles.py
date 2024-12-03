"""
Classifies: CHEBI:48901 thiazoles
"""
from rdkit import Chem

def is_thiazoles(smiles: str):
    """
    Determines if a molecule is a thiazole.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiazole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 5-membered ring
    if not any(len(ring) == 5 for ring in rings.AtomRings()):
        return False, "No 5-membered rings found"

    # Find all aromatic 5-membered rings
    thiazole_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                thiazole_rings.append(ring)

    if not thiazole_rings:
        return False, "No aromatic 5-membered rings found"

    # Check for N and S atoms in the ring
    for ring in thiazole_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        has_nitrogen = any(atom.GetSymbol() == 'N' for atom in atoms)
        has_sulfur = any(atom.GetSymbol() == 'S' for atom in atoms)
        if has_nitrogen and has_sulfur:
            return True, "Thiazole structure found"

    return False, "No thiazole structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48901',
                          'name': 'thiazoles',
                          'definition': 'An azole in which the five-membered '
                                        'heterocyclic aromatic skeleton '
                                        'contains a N atom and one S atom.',
                          'parents': ['CHEBI:38106', 'CHEBI:68452']},
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
    'num_true_positives': 40,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 3,
    'precision': 1.0,
    'recall': 0.9302325581395349,
    'f1': 0.963855421686747,
    'accuracy': None}