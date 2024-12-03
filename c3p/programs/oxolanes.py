"""
Classifies: CHEBI:26912 oxolanes
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_oxolanes(smiles: str):
    """
    Determines if a molecule is an oxolane (tetrahydrofuran or substituted tetrahydrofuran).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxolane, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 5-membered ring
    if not any(len(ring) == 5 for ring in rings.AtomRings()):
        return False, "No 5-membered rings found"

    # Find all 5-membered rings containing an oxygen atom
    oxolane_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in atoms):
                oxolane_rings.append(ring)

    if not oxolane_rings:
        return False, "No 5-membered rings containing an oxygen atom found"

    # Check if the ring is fully saturated (tetrahydrofuran)
    for ring in oxolane_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(atom.GetSymbol() == 'C' or atom.GetSymbol() == 'O' for atom in atoms):
            if all(atom.GetDegree() == 2 for atom in atoms if atom.GetSymbol() == 'C'):
                return True, "Molecule is an oxolane (tetrahydrofuran or substituted tetrahydrofuran)"

    return False, "No fully saturated 5-membered rings containing an oxygen atom found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26912',
                          'name': 'oxolanes',
                          'definition': 'Any oxacycle having an oxolane '
                                        '(tetrahydrofuran) skeleton.',
                          'parents': ['CHEBI:25693', 'CHEBI:38104']},
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
    'num_false_negatives': 90,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}