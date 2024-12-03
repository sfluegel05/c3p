"""
Classifies: CHEBI:24129 furans
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_furans(smiles: str):
    """
    Determines if a molecule is a furan (contains at least one furan ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a furan ring, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the furan substructure
    furan_smiles = "C1=COC=C1"
    furan_mol = Chem.MolFromSmiles(furan_smiles)

    # Check if the molecule contains the furan substructure
    if mol.HasSubstructMatch(furan_mol):
        return True, "Contains at least one furan ring"
    else:
        return False, "Does not contain any furan rings"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24129',
                          'name': 'furans',
                          'definition': 'Compounds containing at least one '
                                        'furan ring.',
                          'parents': ['CHEBI:25693', 'CHEBI:38104']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[20:26:41] SMILES Parse Error: syntax error while parsing: '
             'CC1=CC=C(/C=C\x02/CCN=C2)O1\n'
             '[20:26:41] SMILES Parse Error: Failed parsing SMILES '
             "'CC1=CC=C(/C=C\x02/CCN=C2)O1' for input: "
             "'CC1=CC=C(/C=C\x02/CCN=C2)O1'\n",
    'stdout': '',
    'num_true_positives': 40,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 57,
    'precision': 1.0,
    'recall': 0.41237113402061853,
    'f1': 0.5839416058394161,
    'accuracy': None}