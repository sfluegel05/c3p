"""
Classifies: CHEBI:75769 B vitamin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

# List of SMILES strings for known B vitamins
B_VITAMINS_SMILES = [
    "[H]C(=O)N(C[C@@H]1CNC2=C(N1)C(=O)NC(N)=N2)C1=CC=C(C=C1)C(=O)N[C@@H](CCC(O)=O)C(O)=O",  # (6S)-10-formyltetrahydrofolic acid
    "Nc1nc2N=C[C@H](CNc3ccc(cc3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)N(CO)c2c(=O)[nH]1",  # 5-hydroxymethyl-5,6-dihydrofolic acid
    "[C@H](CCC([O-])=O)(C([O-])=O)NC(C1=CC=C(NCC2N(C3=C(NC2)NC(N)=NC3=O)C=O)C=C1)=O.[Ca+2]",  # Calcium folinate
    "Cc1ncc(CO)c(C([O-])=O)c1O",  # 4-pyridoxate
    "CN1[C@@H](CNC2=CC=C(C=C2)C(=O)N[C@@H](CCC(O)=O)C(O)=O)CNC2=C1C(=O)NC(N)=N2",  # (6S)-5-methyltetrahydrofolic acid
    "[Cl-].CC1=C(CCOP(O)(=O)OP(O)(O)=O)SC=[N+]1CC1=C(N)N=C(C)N=C1",  # thiamine(1+) diphosphate chloride
    "C1(O)=C(C)N=CC(CO)=C1CN.Cl",  # pyridoxamine hydrochloride
    "OC(=O)[C@@H](NC(=O)C1=CC=C(N(CC=2N=C3C(=NC2)N=C(N=C3N)N)C)C=C1)CCC(O)=O.O",  # Methotrexate hydrate
    "C([N+]=1C(=C(SC1)CCO)C)C=2C(=NC(C)=NC2)[NH3+].[Cl-].[Cl-].O.O",  # thiamine hydrochloride dihydrate
    "[Co-3]1234(N5C6=C(C7=[N+]4C(=CC8=[N+]3C(=C(C9=[N+]2[C@@]([C@]5([C@@H]([C@@]6(C)CCC(=O)NC[C@H](OP(O[C@@H]%10[C@H](O[C@H](N%11C=[N+]1C%12=CC(=C(C=C%12%11)C)C)[C@@H]%10O)CO)(=O)[O-])C)CC(=O)N)[H])([C@]([C@@H]9CCC(N)=O)(CC(=O)N)C)C)C)[C@](C)([C@@H]8CCC(=O)N)CC(N)=O)C([C@@H]7CCC(=O)N)(C)C)C)*"  # R-cob(III)alamin
]

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    for b_vitamin_smiles in B_VITAMINS_SMILES:
        b_vitamin_mol = Chem.MolFromSmiles(b_vitamin_smiles)
        if b_vitamin_mol is None:
            continue
        if mol.HasSubstructMatch(b_vitamin_mol):
            return True, "Matches a known B vitamin structure"
    
    return False, "Does not match any known B vitamin structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:75769',
                          'name': 'B vitamin',
                          'definition': 'Any member of the group of eight '
                                        'water-soluble vitamins originally '
                                        'thought to be a single compound '
                                        '(vitamin B) that play important roles '
                                        'in cell metabolism. The group '
                                        'comprises of vitamin B1, B2, B3, B5, '
                                        'B6, B7, B9, and B12 (Around 20 other '
                                        'compounds were once thought to be B '
                                        'vitamins but are no longer classified '
                                        'as such).',
                          'parents': ['CHEBI:35352', 'CHEBI:36963']},
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
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}