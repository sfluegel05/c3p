"""
Classifies: CHEBI:18035 diglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diglyceride, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone
    glycerol_smarts = "[C@@H](O)C(O)CO"
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(glycerol_smarts)):
        return False, "No glycerol backbone found"

    # Count number of ester bonds
    ester_bonds = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O[C@H]"))
    if len(ester_bonds) != 2:
        return False, f"Expected 2 ester bonds, found {len(ester_bonds)}"

    # Check for acyl groups
    acyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)C"))
    if len(acyl_groups) < 2:
        return False, "Less than 2 acyl groups found"

    return True, "Valid diglyceride"

# Example usage:
# print(is_diglyceride("C([C@@](COC(CCCCCCC/C=C\CCCC)=O)(OC(CCCCCCCCCCCCCCCCCCCCCCC)=O)[H])O"))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18035',
                          'name': 'diglyceride',
                          'definition': 'A glyceride that is glycerol in which '
                                        'any two of the hydroxy groups have '
                                        'been acylated. In the structure '
                                        'shown, two of the R groups (positions '
                                        'not specified) are acyl groups while '
                                        'the remaining R group can be either H '
                                        'or an alkyl group.',
                          'parents': ['CHEBI:47778', 'CHEBI:76578']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 95,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}