"""
Classifies: CHEBI:38741 isoflavanones
"""
from rdkit import Chem

def is_isoflavanones(smiles: str):
    """
    Determines if a molecule is an isoflavanone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for isoflavanones
    isoflavanone_smarts = 'O1C(CC=C2)=CC=C2C(=O)C3=C1C=C(O)C=C3'
    
    # Create a molecule from the SMARTS pattern
    isoflavanone_mol = Chem.MolFromSmarts(isoflavanone_smarts)
    
    # Check if the molecule matches the isoflavanone pattern
    if mol.HasSubstructMatch(isoflavanone_mol):
        return True, "Molecule matches isoflavanone pattern"
    
    return False, "Molecule does not match isoflavanone pattern"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38741',
                          'name': 'isoflavanones',
                          'definition': 'Members of the class of isoflavans '
                                        'that have a '
                                        '3,4-dihydro-3-aryl-2H-1-benzopyran-4-one '
                                        'skeleton and its substituted '
                                        'derivatives.',
                          'parents': ['CHEBI:3992', 'CHEBI:72572']},
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
    'num_true_negatives': 11,
    'num_false_negatives': 11,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}