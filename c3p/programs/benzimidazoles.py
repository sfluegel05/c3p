"""
Classifies: CHEBI:22715 benzimidazoles
"""
from rdkit import Chem

def is_benzimidazoles(smiles: str):
    """
    Determines if a molecule is a benzimidazole (benzene ring fused to an imidazole ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzimidazole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring and one 5-membered ring
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]
    five_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 5]

    if not six_membered_rings:
        return False, "No 6-membered rings found"
    if not five_membered_rings:
        return False, "No 5-membered rings found"

    # Check for fused rings
    for six_ring in six_membered_rings:
        for five_ring in five_membered_rings:
            fused_atoms = set(six_ring).intersection(set(five_ring))
            if len(fused_atoms) >= 2:
                # Check if the 5-membered ring is an imidazole
                atoms = [mol.GetAtomWithIdx(i) for i in five_ring]
                if sum(1 for atom in atoms if atom.GetSymbol() == 'N') == 2 and sum(1 for atom in atoms if atom.GetSymbol() == 'C') == 3:
                    return True, "Benzimidazole structure found"
    
    return False, "No benzimidazole structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22715',
                          'name': 'benzimidazoles',
                          'definition': 'An organic heterocyclic compound '
                                        'containing a benzene ring fused to an '
                                        'imidazole ring.',
                          'parents': ['CHEBI:27171', 'CHEBI:38101']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 40,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 0,
    'precision': 0.975609756097561,
    'recall': 1.0,
    'f1': 0.9876543209876543,
    'accuracy': None}