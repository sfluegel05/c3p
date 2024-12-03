"""
Classifies: CHEBI:33299 alkaline earth molecular entity
"""
from rdkit import Chem

def is_alkaline_earth_molecular_entity(smiles: str):
    """
    Determines if a molecule is an alkaline earth molecular entity.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaline earth molecular entity, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Alkaline earth metals: Beryllium (Be), Magnesium (Mg), Calcium (Ca), Strontium (Sr), Barium (Ba), Radium (Ra)
    alkaline_earth_metals = {"Be", "Mg", "Ca", "Sr", "Ba", "Ra"}

    for atom in mol.GetAtoms():
        if atom.GetSymbol() in alkaline_earth_metals:
            return True, f"Contains alkaline earth metal: {atom.GetSymbol()}"

    return False, "No alkaline earth metals found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33299',
                          'name': 'alkaline earth molecular entity',
                          'definition': 'An alkaline earth molecular entity is '
                                        'a molecular entity containing one or '
                                        'more atoms of an alkaline earth '
                                        'metal.',
                          'parents': ['CHEBI:33674']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[23:00:29] SMILES Parse Error: syntax error while parsing: '
             'CCC1=C(C)/C2=C/c3c(C=C)c(C)c4\\C=C5/N=C(C(\\C=C\\C(O)=O)=C/5C)C5=c6c(C(=O)C5C(=O)OC)c(C)c(=CC1=N\x02)n6[Mg]n34\n'
             '[23:00:29] SMILES Parse Error: Failed parsing SMILES '
             "'CCC1=C(C)/C2=C/c3c(C=C)c(C)c4\\C=C5/N=C(C(\\C=C\\C(O)=O)=C/5C)C5=c6c(C(=O)C5C(=O)OC)c(C)c(=CC1=N\x02)n6[Mg]n34' "
             'for input: '
             "'CCC1=C(C)/C2=C/c3c(C=C)c(C)c4\\C=C5/N=C(C(\\C=C\\C(O)=O)=C/5C)C5=c6c(C(=O)C5C(=O)OC)c(C)c(=CC1=N\x02)n6[Mg]n34'\n",
    'stdout': '',
    'num_true_positives': 15,
    'num_false_positives': 0,
    'num_true_negatives': 16,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9375,
    'f1': 0.967741935483871,
    'accuracy': None}