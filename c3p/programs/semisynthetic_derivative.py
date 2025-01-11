"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    A semisynthetic derivative is defined as any organic molecular entity derived from a natural product by partial chemical synthesis.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for common synthetic modifications
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")  # Ester group
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")  # Amide group
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")          # Ether group
    alkyl_halide_pattern = Chem.MolFromSmarts("[Cl,Br,I][CX4]")  # Alkyl halide group

    # Count the number of synthetic modifications
    synthetic_modifications = 0
    if mol.HasSubstructMatch(ester_pattern):
        synthetic_modifications += 1
    if mol.HasSubstructMatch(amide_pattern):
        synthetic_modifications += 1
    if mol.HasSubstructMatch(ether_pattern):
        synthetic_modifications += 1
    if mol.HasSubstructMatch(alkyl_halide_pattern):
        synthetic_modifications += 1

    # If there are multiple synthetic modifications, it is likely a semisynthetic derivative
    if synthetic_modifications >= 2:
        return True, "Contains multiple synthetic modifications (e.g., esters, amides, ethers, alkyl halides)"

    # If there is at least one synthetic modification, it might be a semisynthetic derivative
    if synthetic_modifications == 1:
        return True, "Contains at least one synthetic modification (e.g., ester, amide, ether, alkyl halide)"

    # If no synthetic modifications are found, it is less likely to be a semisynthetic derivative
    return False, "No significant synthetic modifications found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:72588',
                          'name': 'semisynthetic derivative',
                          'definition': 'Any organic molecular entity derived '
                                        'from a natural product by partial '
                                        'chemical synthesis.',
                          'parents': ['CHEBI:50860'],
                          'xrefs': ['Wikipedia:Semisynthesis'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'Cc1ccc(cc1)N=C=Nc1ccc(C)cc1',
                                     'name': '1,3-di(p-tolyl)carbodiimide',
                                     'reason': 'No significant synthetic '
                                               'modifications found'},
                                 {   'smiles': 'N(=CC1=CC=CC=C1)CCC2=CC=CC=C2',
                                     'name': 'N-benzylidene-2-phenylethanamine',
                                     'reason': 'No significant synthetic '
                                               'modifications found'},
                                 {   'smiles': 'C1=CC(=C(C=C1C=C(C#N)C#N)O)[N+](=O)[O-]',
                                     'name': '2-[(3-hydroxy-4-nitrophenyl)methylidene]propanedinitrile',
                                     'reason': 'No significant synthetic '
                                               'modifications found'},
                                 {   'smiles': 'O=C1C2=C(C(=O)C(=C1C(C)C)C)C3[C@@H](C4=CC=CC=C4)C56C2(C(=O)C(=C(O)C5=O)C(C)C)C6[C@H]3C7=CC=CC=C7',
                                     'name': 'Carbocyclinone-534',
                                     'reason': 'No significant synthetic '
                                               'modifications found'},
                                 {   'smiles': 'O=C1C=CC(=O)c2[nH]cnc12',
                                     'name': 'benzimidazole-4,7-quinone',
                                     'reason': 'No significant synthetic '
                                               'modifications found'},
                                 {   'smiles': 'NC(=[NH2+])NCCC[C@H]([NH2+]CCC([O-])=O)C([O-])=O',
                                     'name': 'N(2)-(2-carboxyethyl)-L-arginine '
                                             'dizwitterion',
                                     'reason': 'No significant synthetic '
                                               'modifications found'},
                                 {   'smiles': 'O=C(C1=C(C)C=NC=C1)C',
                                     'name': '4-Acetyl-3-methylpyridine',
                                     'reason': 'No significant synthetic '
                                               'modifications found'},
                                 {   'smiles': 'CC(C)=CCC\\C(C)=C\\Cc1cc(O)ccc1O',
                                     'name': 'geranylhydroquinone',
                                     'reason': 'No significant synthetic '
                                               'modifications found'},
                                 {   'smiles': 'O=C1C(=NNC2=C1C=CC=C2)C(=O)C',
                                     'name': '3-acetyl-4-hydroxycinnoline',
                                     'reason': 'No significant synthetic '
                                               'modifications found'},
                                 {   'smiles': 'CC(C)=CCCC([C@@]1([C@]2([C@H](CC[C@]2([C@@]3(CC[C@]13[H])C)[H])C)[H])[H])=C',
                                     'name': 'prenylkelsoene',
                                     'reason': 'No significant synthetic '
                                               'modifications found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify a molecule as a semisynthetic derivative, we '
                 'need to identify if it is derived from a natural product by '
                 'partial chemical synthesis. This is a challenging task '
                 'because it requires knowledge of the natural product origin '
                 'and the synthetic modifications. However, we can attempt to '
                 'identify common structural features that are often modified '
                 'in semisynthetic derivatives, such as the presence of '
                 'functional groups that are typically introduced through '
                 'chemical synthesis (e.g., esters, amides, ethers, etc.).\n'
                 '\n'
                 'Here is a Python program that attempts to classify a '
                 'molecule as a semisynthetic derivative based on its SMILES '
                 'string:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 43,
    'num_false_positives': 100,
    'num_true_negatives': 14,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3006993006993007,
    'recall': 1.0,
    'f1': 0.46236559139784944,
    'accuracy': 0.3630573248407643,
    'negative_predictive_value': 1.0}