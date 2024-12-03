"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a sphingoid base or ceramide backbone
    sphingoid_smarts = Chem.MolFromSmarts("C[C@H](O)CC=CCCCC")
    ceramide_smarts = Chem.MolFromSmarts("C(=O)N[C@H](CO[C@@H]1O[C@@H](CO)[C@@H](O)[C@@H](O)[C@H]1O)[C@H](O)CC")

    if not mol.HasSubstructMatch(sphingoid_smarts) and not mol.HasSubstructMatch(ceramide_smarts):
        return False, "No sphingoid base or ceramide backbone found"

    # Check for the presence of a glycosidic linkage to O-1 of the sphingoid
    glycosidic_linkage_smarts = Chem.MolFromSmarts("O[C@@H]1O[C@@H](CO)[C@@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glycosidic_linkage_smarts):
        return False, "No glycosidic linkage to O-1 of the sphingoid found"

    return True, "Molecule is a glycosphingolipid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24402',
                          'name': 'glycosphingolipid',
                          'definition': 'A glycosphingolipid is a glycolipid '
                                        'that is a carbohydrate-containing '
                                        'derivative of a sphingoid or '
                                        'ceramide. It is understood that the '
                                        'carbohydrate residue is attached by a '
                                        'glycosidic linkage to O-1 of the '
                                        'sphingoid.',
                          'parents': ['CHEBI:26739', 'CHEBI:33563']},
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
    'num_true_positives': 25,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 96,
    'precision': 0.9615384615384616,
    'recall': 0.2066115702479339,
    'f1': 0.3401360544217687,
    'accuracy': None}