"""
Classifies: CHEBI:35186 terpene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_terpene(smiles: str):
    """
    Determines if a molecule is a terpene.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a terpene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains carbon and hydrogen only
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in {'C', 'H'}:
            return False, "Molecule contains elements other than carbon and hydrogen"

    # Check if the molecule is derived from isoprene units
    isoprene_smiles = "C=C(C)C=C"
    isoprene_mol = Chem.MolFromSmiles(isoprene_smiles)

    if not mol.HasSubstructMatch(isoprene_mol):
        return False, "Molecule does not contain isoprene units"

    return True, "Molecule is a terpene"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35186',
                          'name': 'terpene',
                          'definition': 'A hydrocarbon of biological origin '
                                        'having carbon skeleton formally '
                                        'derived from isoprene '
                                        '[CH2=C(CH3)CH=CH2].',
                          'parents': ['CHEBI:24632', 'CHEBI:24913']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '[23:22:48] SMILES Parse Error: syntax error while parsing: '
             'C1CC(\\C(\\C[C@H]1O)=C/C=C\x02/CC[C@]3([C@](C2)(CC[C@@]3([C@H](C)CCCC(CO)(C)O)[H])[H])C)=C\n'
             '[23:22:48] SMILES Parse Error: Failed parsing SMILES '
             "'C1CC(\\C(\\C[C@H]1O)=C/C=C\x02/CC[C@]3([C@](C2)(CC[C@@]3([C@H](C)CCCC(CO)(C)O)[H])[H])C)=C' "
             'for input: '
             "'C1CC(\\C(\\C[C@H]1O)=C/C=C\x02/CC[C@]3([C@](C2)(CC[C@@]3([C@H](C)CCCC(CO)(C)O)[H])[H])C)=C'\n",
    'stdout': '',
    'num_true_positives': 7,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 45,
    'precision': 1.0,
    'recall': 0.1346153846153846,
    'f1': 0.23728813559322035,
    'accuracy': None}