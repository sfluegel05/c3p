"""
Classifies: CHEBI:36130 cyclic terpene ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_cyclic_terpene_ketone(smiles: str):
    """
    Determines if a molecule is a cyclic terpene ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic terpene ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a ketone group
    ketone_group = Chem.MolFromSmarts('C(=O)')
    if not mol.HasSubstructMatch(ketone_group):
        return False, "No ketone group found"

    # Check for the presence of a carbocyclic ring
    ring_info = mol.GetRingInfo()
    carbocyclic_rings = [ring for ring in ring_info.AtomRings() if all(mol.GetAtomWithIdx(idx).GetSymbol() == 'C' for idx in ring)]
    if not carbocyclic_rings:
        return False, "No carbocyclic ring found"

    # Check for the presence of a terpene skeleton
    terpene_skeletons = [
        Chem.MolFromSmarts('C1CCCCC1'),  # cyclohexane
        Chem.MolFromSmarts('C1CCCC1'),   # cyclopentane
        Chem.MolFromSmarts('C1CCC(C1)C'),# cyclobutane with methyl
        Chem.MolFromSmarts('C1CC(C1)C')  # cyclopropane with ethyl
    ]
    if not any(mol.HasSubstructMatch(ts) for ts in terpene_skeletons):
        return False, "No terpene skeleton found"

    # Check for multiple rings to ensure it's a terpene-like structure
    if len(carbocyclic_rings) < 2:
        return False, "Not enough carbocyclic rings to form a terpene-like structure"

    return True, "Molecule is a cyclic terpene ketone"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36130',
                          'name': 'cyclic terpene ketone',
                          'definition': 'An alicyclic ketone in which the '
                                        'carbocyclic ring structure forms part '
                                        'of a terpene skeleton.',
                          'parents': ['CHEBI:26872', 'CHEBI:36132']},
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
    'num_true_positives': 14,
    'num_false_positives': 6,
    'num_true_negatives': 13,
    'num_false_negatives': 5,
    'precision': 0.7,
    'recall': 0.7368421052631579,
    'f1': 0.717948717948718,
    'accuracy': None}