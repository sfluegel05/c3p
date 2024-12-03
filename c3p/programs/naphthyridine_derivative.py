"""
Classifies: CHEBI:73539 naphthyridine derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_naphthyridine_derivative(smiles: str):
    """
    Determines if a molecule is a naphthyridine derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a naphthyridine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least two rings
    if len(rings.AtomRings()) < 2:
        return False, "Less than two rings found"

    # Check for presence of nitrogen atoms in the rings
    nitrogen_rings = []
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if sum(atom.GetSymbol() == 'N' for atom in atoms) >= 2:
            nitrogen_rings.append(ring)

    if len(nitrogen_rings) < 1:
        return False, "No rings with at least two nitrogen atoms found"

    # Check if the rings are fused to form a naphthyridine core
    for ring1 in nitrogen_rings:
        for ring2 in rings.AtomRings():
            if ring1 != ring2 and len(set(ring1) & set(ring2)) >= 2:
                return True, "Naphthyridine core detected"

    return False, "Naphthyridine core not detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73539',
                          'name': 'naphthyridine derivative',
                          'definition': 'Any organonitrogen heterocyclic '
                                        'compound that is a derivative of a '
                                        'naphthyridine.',
                          'parents': ['CHEBI:27171', 'CHEBI:38101']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 15,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 5,
    'precision': 0.9375,
    'recall': 0.75,
    'f1': 0.8333333333333334,
    'accuracy': None}