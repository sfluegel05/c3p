"""
Classifies: CHEBI:17984 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA (a thioester that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any carboxylic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the coenzyme A (CoA) moiety
    coa_smarts = "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N)C@@H](O)C(=O)NCCC(=O)NCCS"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Check for the presence of the thioester linkage
    thioester_smarts = "C(=O)SCCNC(=O)"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    return True, "Molecule is an acyl-CoA"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17984',
                          'name': 'acyl-CoA',
                          'definition': 'A thioester that results from the '
                                        'formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any carboxylic acid.',
                          'parents': [   'CHEBI:140325',
                                         'CHEBI:20706',
                                         'CHEBI:231540',
                                         'CHEBI:51277']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': 'Python argument types in\n'
             '    Mol.HasSubstructMatch(Mol, NoneType)\n'
             'did not match C++ signature:\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle query, '
             'RDKit::SubstructMatchParameters params=True)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
             'RDKit::SubstructMatchParameters params)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle query, '
             'bool recursionPossible=True, bool useChirality=False, bool '
             'useQueryQueryMatches=False)\n'
             '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
             'bool recursionPossible=True, bool useChirality=False, bool '
             'useQueryQueryMatches=False)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}