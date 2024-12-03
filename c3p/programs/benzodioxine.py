"""
Classifies: CHEBI:64096 benzodioxine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_benzodioxine(smiles: str):
    """
    Determines if a molecule is a benzodioxine (ortho-fused benzene and dioxine rings).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzodioxine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered benzene ring
    benzene_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetSymbol() == 'C' for atom in atoms):
                benzene_rings.append(ring)

    if not benzene_rings:
        return False, "No benzene rings found"

    # Check for dioxine rings (6-membered rings with two oxygens)
    dioxine_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if sum(1 for atom in atoms if atom.GetSymbol() == 'O') == 2:
                dioxine_rings.append(ring)

    if not dioxine_rings:
        return False, "No dioxine rings found"

    # Check for ortho-fusion (shared bond between benzene and dioxine)
    for benzene_ring in benzene_rings:
        for dioxine_ring in dioxine_rings:
            shared_atoms = set(benzene_ring) & set(dioxine_ring)
            if len(shared_atoms) == 2:
                return True, "Benzodioxine structure found"

    return False, "No ortho-fused benzene and dioxine rings found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64096',
                          'name': 'benzodioxine',
                          'definition': 'Any organic heterobicyclic compound '
                                        'containing ortho-fused benzene and '
                                        'dioxine rings.',
                          'parents': ['CHEBI:27171', 'CHEBI:38104']},
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
    'num_true_positives': 14,
    'num_false_positives': 0,
    'num_true_negatives': 14,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}