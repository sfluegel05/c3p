"""
Classifies: CHEBI:50018 glycosyloxyflavone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_glycosyloxyflavone(smiles: str):
    """
    Determines if a molecule is a glycosyloxyflavone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosyloxyflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavone core (C15H10O2)
    flavone_core = Chem.MolFromSmarts('c1cc2c(cc1)oc(=O)cc2')

    # Check if the molecule contains the flavone core
    if not mol.HasSubstructMatch(flavone_core):
        return False, "No flavone core found"

    # Define glycosyl groups (simple sugars)
    glycosyl_groups = [
        Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O'),  # Glucose
        Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O'),  # Galactose
        Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O'),  # Rhamnose
        Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O'),  # Xylose
        Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O')   # Arabinose
    ]

    # Check if the molecule contains at least one glycosyl group
    contains_glycosyl = any(mol.HasSubstructMatch(glycosyl_group) for glycosyl_group in glycosyl_groups)
    if not contains_glycosyl:
        return False, "No glycosyl group found"

    return True, "Molecule is a glycosyloxyflavone"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50018',
                          'name': 'glycosyloxyflavone',
                          'definition': 'A member of the class of flavones '
                                        'having one or more glycosyl residues '
                                        'attached at unspecified positions.',
                          'parents': ['CHEBI:24043', 'CHEBI:24400']},
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
    'num_false_negatives': 22,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}