"""
Classifies: CHEBI:63427 glycosylglycerol derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_glycosylglycerol_derivative(smiles: str):
    """
    Determines if a molecule is a glycosylglycerol derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosylglycerol derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the substructure patterns for glycerol and glycosyl
    glycerol = Chem.MolFromSmarts("C(CO)CO")
    glycosyl = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O1)")

    if glycerol is None or glycosyl is None:
        return False, "Error in defining substructure patterns"

    # Check if the molecule contains a glycerol moiety
    if not mol.HasSubstructMatch(glycerol):
        return False, "No glycerol moiety found"

    # Check if the molecule contains a glycosyl moiety
    if not mol.HasSubstructMatch(glycosyl):
        return False, "No glycosyl moiety found"

    # Check if the glycosyl moiety is attached to the glycerol moiety
    glycosyl_matches = mol.GetSubstructMatches(glycosyl)
    glycerol_matches = mol.GetSubstructMatches(glycerol)

    for glycosyl_match in glycosyl_matches:
        for glycerol_match in glycerol_matches:
            if any(atom in glycosyl_match for atom in glycerol_match):
                return True, "Glycosyl moiety attached to glycerol moiety found"

    return False, "Glycosyl moiety not attached to glycerol moiety"

# Example usage:
# smiles = "CC\C=C/CC1C(CCCCCCCC(=O)OC[C@H](CO[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)OC(=O)CCCCCC2C=CC(=O)C2C\C=C/CC)C=CC1=O"
# result, reason = is_glycosylglycerol_derivative(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63427',
                          'name': 'glycosylglycerol derivative',
                          'definition': 'A glycosyl alditol derivative that is '
                                        'formally obtained from a glycosyl '
                                        'glycerol.',
                          'parents': ['CHEBI:63424']},
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
    'num_true_positives': 21,
    'num_false_positives': 6,
    'num_true_negatives': 14,
    'num_false_negatives': 0,
    'precision': 0.7777777777777778,
    'recall': 1.0,
    'f1': 0.8750000000000001,
    'accuracy': None}