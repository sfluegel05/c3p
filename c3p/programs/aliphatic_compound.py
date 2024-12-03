"""
Classifies: CHEBI:33653 aliphatic compound
"""
from rdkit import Chem

def is_aliphatic_compound(smiles: str):
    """
    Determines if a molecule is an aliphatic compound (any acyclic or cyclic, saturated or unsaturated carbon compound, excluding aromatic compounds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is aromatic
    if mol.GetAromaticAtoms():
        return False, "Molecule contains aromatic atoms"

    # Check if the molecule contains carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
    if not carbon_atoms:
        return False, "Molecule does not contain carbon atoms"

    # Check for the presence of aromatic rings
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            return False, "Molecule contains aromatic rings"

    return True, "Molecule is an aliphatic compound"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33653',
                          'name': 'aliphatic compound',
                          'definition': 'Any acyclic or cyclic, saturated or '
                                        'unsaturated carbon compound, '
                                        'excluding aromatic compounds.',
                          'parents': ['CHEBI:50860']},
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
    'num_true_positives': 44,
    'num_false_positives': 10,
    'num_true_negatives': 10,
    'num_false_negatives': 1,
    'precision': 0.8148148148148148,
    'recall': 0.9777777777777777,
    'f1': 0.8888888888888888,
    'accuracy': None}