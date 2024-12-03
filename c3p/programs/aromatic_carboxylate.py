"""
Classifies: CHEBI:91007 aromatic carboxylate
"""
from rdkit import Chem

def is_aromatic_carboxylate(smiles: str):
    """
    Determines if a molecule is an aromatic carboxylate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic carboxylate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate group [O-]C(=O)
    carboxylate_pattern = Chem.MolFromSmarts('[O-]C(=O)')
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Check for aromatic ring
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]
    if not aromatic_rings:
        return False, "No aromatic rings found"

    # Ensure the carboxylate group is connected to the aromatic ring
    for ring in aromatic_rings:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if mol.GetSubstructMatch(carboxylate_pattern) and neighbor.GetIdx() in mol.GetSubstructMatch(carboxylate_pattern):
                    return True, "Aromatic carboxylate found"
    
    return False, "Carboxylate group not connected to aromatic ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:91007',
                          'name': 'aromatic carboxylate',
                          'definition': 'A carboxylic acic anion obtained by '
                                        'deprotonation of the carboxy group of '
                                        'any aromatic carboxylic acid. Major '
                                        'species at pH 7.3.',
                          'parents': ['CHEBI:29067']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 15,
    'num_false_positives': 0,
    'num_true_negatives': 16,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9375,
    'f1': 0.967741935483871,
    'accuracy': None}