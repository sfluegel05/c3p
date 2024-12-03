"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 5beta-cholanic acid core structure
    cholanic_acid_smarts = Chem.MolFromSmarts("C1[C@@]2([C@]3([C@@]([C@](CC3)([C@@H](CCC(O)=O)C)[H])(CC[C@@]2([C@@]4([C@](C1)(CC(=O)CC4)[H])C)[H])C)[H])[H]")
    if not mol.HasSubstructMatch(cholanic_acid_smarts):
        return False, "Does not contain 5beta-cholanic acid core structure"

    # Check for hydroxy groups
    hydroxy_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]"))
    if len(hydroxy_groups) == 0:
        return False, "No hydroxy groups found"

    # Check for amides with glycine or taurine
    glycine_amide_smarts = Chem.MolFromSmarts("NCC(=O)[O-]")
    taurine_amide_smarts = Chem.MolFromSmarts("NCCS(=O)(=O)[O-]")
    if not mol.HasSubstructMatch(glycine_amide_smarts) and not mol.HasSubstructMatch(taurine_amide_smarts):
        return False, "No glycine or taurine amide found"

    return True, "Molecule is a bile acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:3098',
                          'name': 'bile acid',
                          'definition': 'Any member of a group of '
                                        'hydroxy-5beta-cholanic acids occuring '
                                        'in bile, where they are present as '
                                        'the sodium salts of their amides with '
                                        'glycine or taurine. In mammals bile '
                                        'acids almost invariably have '
                                        '5beta-configuration.',
                          'parents': ['CHEBI:138366', 'CHEBI:24663']},
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
    'num_true_negatives': 20,
    'num_false_negatives': 42,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}