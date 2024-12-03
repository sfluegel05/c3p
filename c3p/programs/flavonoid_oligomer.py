"""
Classifies: CHEBI:72720 flavonoid oligomer
"""
from rdkit import Chem

def is_flavonoid_oligomer(smiles: str):
    """
    Determines if a molecule is a flavonoid oligomer.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid oligomer, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least two 6-membered rings
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings found"

    # Find all aromatic 6-membered rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic 6-membered rings found"

    # Check if the molecule contains benzopyran units
    benzopyran_pattern = Chem.MolFromSmarts("c1cc2occc2cc1")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No benzopyran units found"

    # Check if the molecule has multiple benzopyran units
    benzopyran_matches = mol.GetSubstructMatches(benzopyran_pattern)
    if len(benzopyran_matches) < 2:
        return False, "Less than two benzopyran units found"

    # Check for the presence of phenylpropanoid units
    phenylpropanoid_pattern = Chem.MolFromSmarts("c1ccccc1CCC")
    if not mol.HasSubstructMatch(phenylpropanoid_pattern):
        return False, "No phenylpropanoid units found"

    return True, "Molecule is a flavonoid oligomer"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:72720',
                          'name': 'flavonoid oligomer',
                          'definition': 'A phenylpropanoid obtained by the '
                                        'coupling of two or more units of '
                                        'aryl-substituted benzopyran units and '
                                        'their substituted derivatives.',
                          'parents': [   'CHEBI:26004',
                                         'CHEBI:38443',
                                         'CHEBI:72544']},
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
    'num_true_negatives': 14,
    'num_false_negatives': 14,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}