"""
Classifies: CHEBI:35348 3beta-sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_3beta_sterol(smiles: str):
    """
    Determines if a molecule is a 3beta-sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-sterol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a sterol scaffold
    sterol_substructure = Chem.MolFromSmarts('C1CCC2C3CCC4CC(O)CCC4C3CCC12')
    if not mol.HasSubstructMatch(sterol_substructure):
        return False, "No sterol scaffold found"

    # Check for the 3beta-hydroxy group
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    for center in chiral_centers:
        idx, chirality = center
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == 'C' and chirality == 'S':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(idx, neighbor.GetIdx()).GetBondType().name == 'SINGLE':
                    return True, "3beta-sterol with 3beta-hydroxy group"

    return False, "No 3beta-hydroxy group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35348',
                          'name': '3beta-sterol',
                          'definition': 'A sterol in which the hydroxy group '
                                        'at position 3 has beta- '
                                        'configuration.',
                          'parents': ['CHEBI:15889', 'CHEBI:36836']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 3,
    'num_false_positives': 1,
    'num_true_negatives': 13,
    'num_false_negatives': 11,
    'precision': 0.75,
    'recall': 0.21428571428571427,
    'f1': 0.3333333333333333,
    'accuracy': None}