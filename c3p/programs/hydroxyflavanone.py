"""
Classifies: CHEBI:24697 hydroxyflavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_hydroxyflavanone(smiles: str):
    """
    Determines if a molecule is a hydroxyflavanone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for flavanone core structure
    flavanone_smarts = 'O=C1CC[C@H](Oc2c1ccccc2)c1ccccc1'
    flavanone = Chem.MolFromSmarts(flavanone_smarts)
    if not mol.HasSubstructMatch(flavanone):
        return False, "No flavanone core structure found"

    # Check for at least one hydroxy group
    hydroxy_smarts = '[OH]'
    hydroxy = Chem.MolFromSmarts(hydroxy_smarts)
    if not mol.HasSubstructMatch(hydroxy):
        return False, "No hydroxy groups found"

    # Check if hydroxy groups are attached to the flavanone core
    hydroxy_matches = mol.GetSubstructMatches(hydroxy)
    flavanone_matches = mol.GetSubstructMatches(flavanone)
    flavanone_atoms = set(atom for match in flavanone_matches for atom in match)
    
    for match in hydroxy_matches:
        if any(atom in flavanone_atoms for atom in match):
            return True, "Molecule is a hydroxyflavanone"

    return False, "Hydroxy groups are not attached to the flavanone core"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24697',
                          'name': 'hydroxyflavanone',
                          'definition': 'A member of the class of flavanones '
                                        'that consists of flavanone with one '
                                        'or more hydroxy substituents.',
                          'parents': ['CHEBI:28863', 'CHEBI:33822']},
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