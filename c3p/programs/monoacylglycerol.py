"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for acyl group
    acyl_pattern = Chem.MolFromSmarts("C(=O)[#6]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "No acyl group found"

    # Check that only one acyl group is attached to the glycerol backbone
    acyl_count = 0
    for match in acyl_matches:
        acyl_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in acyl_atom.GetNeighbors():
            if neighbor.GetSmarts() == "O" and neighbor.GetNeighbors()[0].GetSmarts() == "C(O)CO":
                acyl_count += 1

    if acyl_count != 1:
        return False, "More than one acyl group found"

    return True, "Valid monoacylglycerol"

# Example usage
smiles = "CCCCCCCCCCCC(=O)OCC(O)CO"
print(is_monoacylglycerol(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17408',
                          'name': 'monoacylglycerol',
                          'definition': 'A glyceride in which any one of the R '
                                        'groups (position not specified) is an '
                                        'acyl group while the remaining two R '
                                        'groups can be either H or alkyl '
                                        'groups.',
                          'parents': ['CHEBI:47778']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(False, 'More than one acyl group found')\n",
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 18,
    'num_false_negatives': 18,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}