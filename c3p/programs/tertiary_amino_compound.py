"""
Classifies: CHEBI:50996 tertiary amino compound
"""
from rdkit import Chem

def is_tertiary_amino_compound(smiles: str):
    """
    Determines if a molecule is a tertiary amino compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amino compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetDegree() == 3:
            # Check if the nitrogen atom is bonded to three carbon atoms
            if all(neigh.GetSymbol() == 'C' for neigh in atom.GetNeighbors()):
                return True, "Tertiary amino compound found"
    
    return False, "No tertiary amino compound found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50996',
                          'name': 'tertiary amino compound',
                          'definition': 'A compound formally derived from '
                                        'ammonia by replacing three hydrogen '
                                        'atoms by organyl groups.',
                          'parents': ['CHEBI:50047']},
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
    'num_true_positives': 138,
    'num_false_positives': 7,
    'num_true_negatives': 13,
    'num_false_negatives': 2,
    'precision': 0.9517241379310345,
    'recall': 0.9857142857142858,
    'f1': 0.968421052631579,
    'accuracy': None}