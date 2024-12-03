"""
Classifies: CHEBI:64898 anionic amino-acid residue
"""
from rdkit import Chem

def is_anionic_amino_acid_residue(smiles: str):
    """
    Determines if a molecule is an anionic amino-acid residue.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anionic amino-acid residue, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of amino acid residue backbone (N-C-C(=O)-)
    amino_acid_residue = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()):
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    for sub_neighbor in neighbor.GetNeighbors():
                        if sub_neighbor.GetSymbol() == 'C' and any(sub_sub_neighbor.GetSymbol() == 'O' for sub_sub_neighbor in sub_neighbor.GetNeighbors()):
                            amino_acid_residue = True
                            break
                if amino_acid_residue:
                    break
        if amino_acid_residue:
            break

    if not amino_acid_residue:
        return False, "No amino-acid residue backbone found"

    # Check for overall negative charge
    charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if charge >= 0:
        return False, "Molecule does not have an overall negative charge"

    return True, "Molecule is an anionic amino-acid residue"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64898',
                          'name': 'anionic amino-acid residue',
                          'definition': 'An amino-acid residue carrying an '
                                        'overall negative charge.',
                          'parents': ['CHEBI:64775']},
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
    'num_true_positives': 39,
    'num_false_positives': 9,
    'num_true_negatives': 11,
    'num_false_negatives': 3,
    'precision': 0.8125,
    'recall': 0.9285714285714286,
    'f1': 0.8666666666666666,
    'accuracy': None}