"""
Classifies: CHEBI:55370 imidazolidinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_imidazolidinone(smiles: str):
    """
    Determines if a molecule is an imidazolidinone (An imidazolidine containing one or more oxo groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an imidazolidinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for imidazolidinone
    imidazolidinone_pattern = Chem.MolFromSmarts('C1NC(=O)NC1')

    if not mol.HasSubstructMatch(imidazolidinone_pattern):
        return False, "No imidazolidinone core structure found"

    # Check for oxo groups
    oxo_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetTotalDegree() == 1]

    if not oxo_groups:
        return False, "No oxo groups found"

    return True, "Valid imidazolidinone"

# Examples for testing
smiles_list = [
    "C(N1C(=O)N(CC1=O)/N=C/C=2OC(=CC2)[N+](=O)[O-])N3CCOCC3",  # nifurfoline
    "C(CN1C(NC([C@H]1CCCCCCC(O)=O)=O)=O)[C@@H](C2CCCCC2)O",  # (3S,4R)-BW 245C
    "O=C1N(N=O)CCN1N=O",  # 1,3-Dinitroso-2-imidazolidinone
    "CCN1C(=O)C(=CC2=CC(=C(C(=C2)OCC)OCC3=CC=C(C=C3)C(=O)O)CC=C)NC1=O",  # 4-[[2-ethoxy-4-[(1-ethyl-2,5-dioxo-4-imidazolidinylidene)methyl]-6-prop-2-enylphenoxy]methyl]benzoic acid
    "O=C1N(C(=O)NC1C2=CC=C(O)C=C2)CC",  # p-Hydroxyl-ethotoin
    "CC1(C(=O)N(C(=O)N1)CC(=O)N2CCOCC2)C3=CC=CC=C3",  # 5-methyl-3-[2-(4-morpholinyl)-2-oxoethyl]-5-phenylimidazolidine-2,4-dione
    "OC(=O)C1NC(=O)NC1",  # 2-Imidazolidone-4-carboxylic acid
    "O=C1NC(=O)C(=O)N1CCc1ccccc1",  # 1-(2-phenylethyl)imidazolidine-2,4,5-trione
    "O=C1N(CCN1C)C",  # 1,3-Dimethyl-2-imidazolidinon
    "CC(C)=CC1C(C(=O)OCN2C(=O)CN(CC#C)C2=O)C1(C)C",  # imiprothrin
    "O=C1N(C(=O)NC1(CC)C2=CC=C(O)C=C2)C",  # Hydroxymephenytoin
    "O=C1NC(=O)N[C@@]1(CC)C2=CC=CC=C2"  # S-nirvanol
]

for smiles in smiles_list:
    result, reason = is_imidazolidinone(smiles)
    print(f"SMILES: {smiles} -> {result}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:55370',
                          'name': 'imidazolidinone',
                          'definition': 'An imidazolidine containing one or '
                                        'more oxo groups.',
                          'parents': ['CHEBI:38261']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: C(N1C(=O)N(CC1=O)/N=C/C=2OC(=CC2)[N+](=O)[O-])N3CCOCC3 '
              '-> True, Valid imidazolidinone\n'
              'SMILES: C(CN1C(NC([C@H]1CCCCCCC(O)=O)=O)=O)[C@@H](C2CCCCC2)O -> '
              'True, Valid imidazolidinone\n'
              'SMILES: O=C1N(N=O)CCN1N=O -> True, Valid imidazolidinone\n'
              'SMILES: '
              'CCN1C(=O)C(=CC2=CC(=C(C(=C2)OCC)OCC3=CC=C(C=C3)C(=O)O)CC=C)NC1=O '
              '-> True, Valid imidazolidinone\n'
              'SMILES: O=C1N(C(=O)NC1C2=CC=C(O)C=C2)CC -> True, Valid '
              'imidazolidinone\n'
              'SMILES: CC1(C(=O)N(C(=O)N1)CC(=O)N2CCOCC2)C3=CC=CC=C3 -> True, '
              'Valid imidazolidinone\n'
              'SMILES: OC(=O)C1NC(=O)NC1 -> True, Valid imidazolidinone\n'
              'SMILES: O=C1NC(=O)C(=O)N1CCc1ccccc1 -> True, Valid '
              'imidazolidinone\n'
              'SMILES: O=C1N(CCN1C)C -> True, Valid imidazolidinone\n'
              'SMILES: CC(C)=CC1C(C(=O)OCN2C(=O)CN(CC#C)C2=O)C1(C)C -> True, '
              'Valid imidazolidinone\n'
              'SMILES: O=C1N(C(=O)NC1(CC)C2=CC=C(O)C=C2)C -> True, Valid '
              'imidazolidinone\n'
              'SMILES: O=C1NC(=O)N[C@@]1(CC)C2=CC=CC=C2 -> True, Valid '
              'imidazolidinone\n',
    'num_true_positives': 12,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}