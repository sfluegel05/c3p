"""
Classifies: CHEBI:59412 amino tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is an amino tetrasaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino tetrasaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has 4 saccharide units
    num_saccharides = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetTotalDegree() == 2:
            # Oxygen atom connected to two carbons, likely part of a glycosidic bond
            num_saccharides += 1

    if num_saccharides < 4:
        return False, "Less than 4 saccharide units detected"

    # Check for the presence of amino groups
    has_amino_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetDegree() > 1:
            has_amino_group = True
            break

    if not has_amino_group:
        return False, "No amino group detected"

    return True, "Molecule is an amino tetrasaccharide"

# Example usage
smiles = "[H][C@@]1(O[C@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2NC(C)=O)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@@H]1O)[C@@H](O)CO"
print(is_amino_tetrasaccharide(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59412',
                          'name': 'amino tetrasaccharide',
                          'definition': 'A tetrasaccharide derivative having '
                                        'one or more substituted or '
                                        'unsubstituted amino groups in place '
                                        'of hydroxy groups at unspecified '
                                        'positions.',
                          'parents': ['CHEBI:22483', 'CHEBI:63567']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Molecule is an amino tetrasaccharide')\n",
    'num_true_positives': 22,
    'num_false_positives': 4,
    'num_true_negatives': 16,
    'num_false_negatives': 0,
    'precision': 0.8461538461538461,
    'recall': 1.0,
    'f1': 0.9166666666666666,
    'accuracy': None}