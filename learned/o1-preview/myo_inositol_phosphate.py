"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: CHEBI:12348 myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component has myo-configuration.
    This involves a cyclohexane ring with six hydroxyl or phosphate groups attached.
    Due to challenges in matching stereochemistry, this implementation focuses on structural features.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find cyclohexane rings
    rings = mol.GetRingInfo().AtomRings()
    cyclohexane_rings = []
    for ring in rings:
        if len(ring) == 6:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            if all(atom.GetAtomicNum() == 6 for atom in atoms_in_ring):
                cyclohexane_rings.append(ring)

    if not cyclohexane_rings:
        return False, "No cyclohexane ring found"

    # For each cyclohexane ring, check for substitution
    for ring in cyclohexane_rings:
        phosphate_found = False
        all_oxygen_substituents = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Get neighbors that are not in the ring
            neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring]
            if not neighbors:
                all_oxygen_substituents = False
                break  # Unsubstituted carbon atom in ring
            else:
                oxygen_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 8]
                if not oxygen_neighbors:
                    all_oxygen_substituents = False
                    break  # Substituent is not oxygen

                for oxy in oxygen_neighbors:
                    # Check if oxygen is connected to a phosphorus atom
                    oxy_neighbors = oxy.GetNeighbors()
                    for oxy_nbr in oxy_neighbors:
                        if oxy_nbr.GetAtomicNum() == 15:
                            phosphate_found = True
                            break  # Phosphate group found
        if all_oxygen_substituents and phosphate_found:
            return True, "Molecule contains cyclohexane ring with oxygen substituents and at least one phosphate group"
    
    return False, "Molecule does not match myo-inositol phosphate structure"
    
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:12348',
        'name': 'myo-inositol phosphate',
        'definition': 'An inositol phosphate in which the inositol component has myo-configuration.',
        'parents': ['CHEBI:24848', 'CHEBI:24845']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}