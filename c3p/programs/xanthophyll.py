"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll (a subclass of carotenoids consisting of the oxygenated carotenes).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is a carotenoid (long conjugated chain with alternating single and double bonds)
    chain_length = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            chain_length += 1

    if chain_length < 9:  # Carotenoids typically have long conjugated chains
        return False, "Not a carotenoid (insufficient conjugated double bonds)"

    # Check for oxygen atoms indicating oxygenated carotenes
    oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O']
    if not oxygen_atoms:
        return False, "No oxygen atoms found (not oxygenated)"

    return True, "Molecule is a xanthophyll (oxygenated carotene)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27325',
                          'name': 'xanthophyll',
                          'definition': 'A subclass of carotenoids consisting '
                                        'of the oxygenated carotenes.',
                          'parents': ['CHEBI:23044']},
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
    'num_true_positives': 22,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 0.9166666666666666,
    'recall': 1.0,
    'f1': 0.9565217391304348,
    'accuracy': None}