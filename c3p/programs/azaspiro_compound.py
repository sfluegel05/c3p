"""
Classifies: CHEBI:35624 azaspiro compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_azaspiro_compound(smiles: str):
    """
    Determines if a molecule is an azaspiro compound.
    An azaspiro compound is a spiro compound in which at least one of the cyclic components is a nitrogen heterocycle.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azaspiro compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check if the molecule has at least one spiro center
    spiro_centers = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    if spiro_centers == 0:
        return False, "No spiro centers found"

    # Check if at least one of the cyclic components is a nitrogen heterocycle
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(atom.GetSymbol() == 'N' for atom in atoms):
            return True, "Contains a nitrogen heterocycle in a spiro compound"

    return False, "No nitrogen heterocycles found in the spiro compound"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35624',
                          'name': 'azaspiro compound',
                          'definition': 'An azaspiro compound is a spiro '
                                        'compound in which at least one of the '
                                        'cyclic components is a nitrogen '
                                        'heterocyle.',
                          'parents': ['CHEBI:33599', 'CHEBI:38101']},
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
    'num_true_positives': 12,
    'num_false_positives': 2,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 0.8571428571428571,
    'recall': 1.0,
    'f1': 0.923076923076923,
    'accuracy': None}