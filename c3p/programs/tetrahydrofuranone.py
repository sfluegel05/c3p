"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone (Any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
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

    # Find all 5-membered rings
    five_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 5]

    for ring in five_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # Check if the ring is an oxolane (tetrahydrofuran)
        if any(atom.GetSymbol() == 'O' for atom in atoms):
            # Check for oxo- substituent (C=O) on the ring
            for atom in atoms:
                if atom.GetSymbol() == 'C':
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            return True, "Tetrahydrofuranone detected"
    return False, "No tetrahydrofuranone structure detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47016',
                          'name': 'tetrahydrofuranone',
                          'definition': 'Any oxolane having an oxo- '
                                        'substituent at any position on the '
                                        'tetrahydrofuran ring.',
                          'parents': ['CHEBI:26912']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 11,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}