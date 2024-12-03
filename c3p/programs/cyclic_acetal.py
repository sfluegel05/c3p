"""
Classifies: CHEBI:59770 cyclic acetal
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_cyclic_acetal(smiles: str):
    """
    Determines if a molecule is a cyclic acetal.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic acetal, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for rings with at least one acetal carbon (a carbon bonded to two oxygens)
    for ring in rings.AtomRings():
        oxygen_count = 0
        acetal_carbon_found = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                # Check if this carbon is bonded to two oxygens
                oxygen_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']
                if len(oxygen_neighbors) == 2:
                    acetal_carbon_found = True
                    oxygen_count += 2
            elif atom.GetSymbol() == 'O':
                oxygen_count += 1

        if acetal_carbon_found and oxygen_count >= 2:
            return True, "Cyclic acetal found"

    return False, "No cyclic acetal found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59770',
                          'name': 'cyclic acetal',
                          'definition': 'An acetal in the molecule of which '
                                        'the acetal carbon and one or both '
                                        'oxygen atoms thereon are members of a '
                                        'ring.',
                          'parents': ['CHEBI:24532', 'CHEBI:59769']},
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
    'num_true_positives': 9,
    'num_false_positives': 6,
    'num_true_negatives': 4,
    'num_false_negatives': 1,
    'precision': 0.6,
    'recall': 0.9,
    'f1': 0.7200000000000001,
    'accuracy': None}