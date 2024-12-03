"""
Classifies: CHEBI:24922 isoquinolines
"""
from rdkit import Chem

def is_isoquinolines(smiles: str):
    """
    Determines if a molecule is an isoquinoline or its substitution derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoquinoline or its substitution derivatives, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the isoquinoline substructure
    isoquinoline_smarts = "c1ccc2ncccc2c1"
    isoquinoline_substructure = Chem.MolFromSmarts(isoquinoline_smarts)
    
    # Define the 3,4-dihydroisoquinoline substructure
    dihydroisoquinoline_smarts = "C1CNc2ccccc2C1"
    dihydroisoquinoline_substructure = Chem.MolFromSmarts(dihydroisoquinoline_smarts)

    if mol.HasSubstructMatch(isoquinoline_substructure) or mol.HasSubstructMatch(dihydroisoquinoline_substructure):
        return True, "Isoquinoline or its substitution derivatives found"
    else:
        return False, "No isoquinoline core found"

# Example usage:
# smiles = "COc1cc(C[C@@H]2NCCc3cc(OC)c(O)cc23)ccc1O"
# result, reason = is_isoquinolines(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24922',
                          'name': 'isoquinolines',
                          'definition': 'A class of organic heteropolycyclic '
                                        'compound consisting of isoquinoline '
                                        'and its substitution derivatives.',
                          'parents': ['CHEBI:38101', 'CHEBI:38166']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 3,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 63,
    'precision': 1.0,
    'recall': 0.045454545454545456,
    'f1': 0.08695652173913045,
    'accuracy': None}