"""
Classifies: CHEBI:140326 tertiary carboxamide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_carboxamide(smiles: str):
    """
    Determines if a molecule is a tertiary carboxamide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary carboxamide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for tertiary carboxamide
    tertiary_carboxamide_smarts = "[NX3][CX3](=O)[NX3]([CX4])[CX4]"
    
    # Create a molecule pattern object
    pattern = Chem.MolFromSmarts(tertiary_carboxamide_smarts)
    
    if pattern is None:
        return None, None

    # Check if the molecule matches the tertiary carboxamide pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains a tertiary carboxamide group"
    else:
        return False, "Molecule does not contain a tertiary carboxamide group"

# Example usage:
# smiles = "CC1=CC=CC=C1C(=O)NC1=CC(C)=C(C=C1)C(=O)N1CCC[C@@H](O)C2=CC(Cl)=CC=C12"
# result, reason = is_tertiary_carboxamide(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:140326',
                          'name': 'tertiary carboxamide',
                          'definition': 'A carboxamide resulting from the '
                                        'formal condensation of a carboxylic '
                                        'acid with a secondary amine; formula '
                                        'RC(=O)NHR(1)R(2).',
                          'parents': ['CHEBI:37622']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 59-60: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}