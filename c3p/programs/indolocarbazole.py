"""
Classifies: CHEBI:51915 indolocarbazole
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_indolocarbazole(smiles: str):
    """
    Determines if a molecule is an indolocarbazole (based upon an indolo[2,3-a]carbazole skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indolocarbazole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the indolo[2,3-a]carbazole core substructure
    indolocarbazole_smiles = "C1=CC=C2C(=C1)C3=NC4=CC=CC=C4C=C3N2"
    indolocarbazole_mol = Chem.MolFromSmiles(indolocarbazole_smiles)
    
    if indolocarbazole_mol is None:
        return None, None  # Shouldn't happen, but just in case the core structure is invalid

    # Check if the molecule contains the indolocarbazole core
    if mol.HasSubstructMatch(indolocarbazole_mol):
        return True, "Molecule contains the indolo[2,3-a]carbazole core structure"
    else:
        return False, "Molecule does not contain the indolo[2,3-a]carbazole core structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51915',
                          'name': 'indolocarbazole',
                          'definition': 'Compounds based upon an '
                                        'indolo[2,3-a]carbazole skeleton.',
                          'parents': ['CHEBI:38166']},
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
    'num_true_negatives': 13,
    'num_false_negatives': 13,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}