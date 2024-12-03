"""
Classifies: CHEBI:23677 diazole
"""
from rdkit import Chem

def is_diazole(smiles: str):
    """
    Determines if a molecule is a diazole (an azole that is either one of a pair of heterocyclic organic compounds comprising three carbon atoms and two nitrogen atoms arranged in a ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diazole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for 5-membered rings
    five_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 5]

    if not five_membered_rings:
        return False, "No 5-membered rings found"

    # Check for rings with 3 carbon atoms and 2 nitrogen atoms
    for ring in five_membered_rings:
        carbon_count = 0
        nitrogen_count = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                carbon_count += 1
            elif atom.GetSymbol() == 'N':
                nitrogen_count += 1
        if carbon_count == 3 and nitrogen_count == 2:
            return True, "Molecule is a diazole"

    return False, "No rings with 3 carbon atoms and 2 nitrogen atoms found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23677',
                          'name': 'diazole',
                          'definition': 'An azole that is either one of a pair '
                                        'of heterocyclic organic compounds '
                                        'comprising three carbon atoms and two '
                                        'nitrogen atoms arranged in a ring.',
                          'parents': ['CHEBI:68452']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 24-25: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}