"""
Classifies: CHEBI:140324 primary carboxamide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_primary_carboxamide(smiles: str):
    """
    Determines if a molecule is a primary carboxamide (a carboxamide resulting from the formal condensation of a carboxylic acid with ammonia).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary carboxamide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the primary carboxamide group
    carboxamide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(carboxamide_pattern):
        return False, "No carboxamide group (C(=O)N) found"

    # Check if the nitrogen in the carboxamide group has exactly two hydrogens (i.e., primary amide)
    for match in mol.GetSubstructMatches(carboxamide_pattern):
        carbonyl_carbon = mol.GetAtomWithIdx(match[0])
        amide_nitrogen = mol.GetAtomWithIdx(match[1])

        # Check if the carbon is double-bonded to oxygen
        if not any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetOtherAtom(carbonyl_carbon).GetSymbol() == 'O' for bond in carbonyl_carbon.GetBonds()):
            continue

        # Check if the nitrogen is bonded to exactly two hydrogens
        hydrogen_count = sum(1 for neighbor in amide_nitrogen.GetNeighbors() if neighbor.GetSymbol() == 'H')
        if hydrogen_count == 2:
            return True, "Primary carboxamide group found"

    return False, "No primary carboxamide group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:140324',
                          'name': 'primary carboxamide',
                          'definition': 'A carboxamide resulting from the '
                                        'formal condensation of a carboxylic '
                                        'acid with ammonia; formula RC(=O)NH2.',
                          'parents': ['CHEBI:37622']},
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
    'num_true_negatives': 10,
    'num_false_negatives': 10,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}