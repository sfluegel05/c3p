"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered benzopyran ring
    benzopyran_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in atoms):
                benzopyran_rings.append(ring)

    if not benzopyran_rings:
        return False, "No benzopyran rings found"

    # Check for aryl substituent at position 2 of the benzopyran ring
    for ring in benzopyran_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        oxygen_idx = next(i for i, atom in enumerate(atoms) if atom.GetSymbol() == 'O')
        position_2_idx = (oxygen_idx + 1) % len(atoms)
        position_2_atom = atoms[position_2_idx]

        if not any(neighbor.GetSymbol() == 'C' and neighbor.GetIsAromatic() for neighbor in position_2_atom.GetNeighbors()):
            return False, "No aryl substituent at position 2 of the benzopyran ring"

    return True, "Molecule is a flavonoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47916',
                          'name': 'flavonoid',
                          'definition': "Any member of the 'superclass' "
                                        'flavonoids whose skeleton is based on '
                                        '1-benzopyran with an aryl substituent '
                                        'at position 2. The term was '
                                        'originally restricted to natural '
                                        'products, but is now also used to '
                                        'describe semi-synthetic and fully '
                                        'synthetic compounds.',
                          'parents': ['CHEBI:72544']},
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
    'num_true_positives': 109,
    'num_false_positives': 11,
    'num_true_negatives': 9,
    'num_false_negatives': 67,
    'precision': 0.9083333333333333,
    'recall': 0.6193181818181818,
    'f1': 0.7364864864864865,
    'accuracy': None}