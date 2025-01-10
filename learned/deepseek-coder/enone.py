"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: CHEBI:51689 enone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha,beta-unsaturated ketone with the general formula R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) =/= H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the enone pattern: C=C-C=O where the carbonyl carbon is not bonded to hydrogen
    enone_pattern = Chem.MolFromSmarts("[CX3]=[CX3]-[CX3]=[OX1]")
    if not mol.HasSubstructMatch(enone_pattern):
        return False, "No alpha,beta-unsaturated ketone pattern found"

    # Verify that the carbonyl carbon is not bonded to hydrogen
    for match in mol.GetSubstructMatches(enone_pattern):
        carbonyl_carbon_idx = match[2]  # Index of the carbonyl carbon in the pattern
        carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_idx)
        if carbonyl_carbon.GetTotalNumHs() > 0:
            return False, "Carbonyl carbon is bonded to hydrogen, not an enone"

    # Check for conjugation between the double bond and the carbonyl group
    # This is already implied by the pattern, but we can double-check
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                # Check if the double bond is conjugated to a carbonyl group
                for neighbor in begin_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalValence() == 2:
                        return True, "Contains alpha,beta-unsaturated ketone (enone) structure"
                for neighbor in end_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalValence() == 2:
                        return True, "Contains alpha,beta-unsaturated ketone (enone) structure"

    return False, "No conjugated alpha,beta-unsaturated ketone structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51689',
                          'name': 'enone',
                          'definition': 'An alpha,beta-unsaturated ketone of '
                                        'general formula R(1)R(2)C=CR(3)-C(=O)R(4) '
                                        '(R(4) =/= H) in which the C=O function '
                                        'is conjugated to a C=C double bond at '
                                        'the alpha,beta position.',
                          'parents': ['CHEBI:51686', 'CHEBI:51687']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}