"""
Classifies: CHEBI:74927 furopyran
"""
from rdkit import Chem

def is_furopyran(smiles: str):
    """
    Determines if a molecule is a furopyran (containing ortho-fused furan and pyran rings).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furopyran, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo().AtomRings()

    # Identify 5-membered and 6-membered rings
    five_membered_rings = [ring for ring in rings if len(ring) == 5]
    six_membered_rings = [ring for ring in rings if len(ring) == 6]

    if not five_membered_rings or not six_membered_rings:
        return False, "No 5-membered or 6-membered rings found"

    # Check for ortho-fused furan and pyran rings
    for five_ring in five_membered_rings:
        for six_ring in six_membered_rings:
            # Check if the rings share two atoms (ortho-fusion)
            shared_atoms = set(five_ring) & set(six_ring)
            if len(shared_atoms) == 2:
                # Check if 5-membered ring is a furan (contains an oxygen)
                if any(mol.GetAtomWithIdx(atom).GetSymbol() == 'O' for atom in five_ring):
                    # Check if 6-membered ring is a pyran (contains an oxygen)
                    if any(mol.GetAtomWithIdx(atom).GetSymbol() == 'O' for atom in six_ring):
                        return True, "Molecule contains ortho-fused furan and pyran rings"
    
    return False, "No ortho-fused furan and pyran rings found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:74927',
                          'name': 'furopyran',
                          'definition': 'Any organic heterobicyclic compound '
                                        'containing ortho-fused furan and '
                                        'pyran rings.',
                          'parents': ['CHEBI:27171', 'CHEBI:38104']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 28-29: malformed \\N character escape (<string>, line '
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