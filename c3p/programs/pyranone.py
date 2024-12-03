"""
Classifies: CHEBI:37963 pyranone
"""
from rdkit import Chem

def is_pyranone(smiles: str):
    """
    Determines if a molecule is a pyranone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyranone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]
    if not six_membered_rings:
        return False, "No 6-membered rings found"

    for ring in six_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        oxygen_count = sum(1 for atom in atoms if atom.GetSymbol() == 'O')
        carbonyl_count = sum(1 for atom in atoms if atom.GetSymbol() == 'C' and any(bond.GetBondTypeAsDouble() == 2.0 and bond.GetOtherAtom(atom).GetSymbol() == 'O' for bond in atom.GetBonds()))

        if oxygen_count == 1 and carbonyl_count == 1:
            return True, "Molecule is a pyranone"

    return False, "No suitable 6-membered ring with one oxygen and one oxo substituent found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37963',
                          'name': 'pyranone',
                          'definition': 'Any of a class of cyclic chemical '
                                        'compounds that contain an unsaturated '
                                        'six-membered ring with one ring '
                                        'oxygen atom and an oxo substituent.',
                          'parents': ['CHEBI:26407']},
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
    'num_true_positives': 78,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 2,
    'precision': 0.975,
    'recall': 0.975,
    'f1': 0.975,
    'accuracy': None}