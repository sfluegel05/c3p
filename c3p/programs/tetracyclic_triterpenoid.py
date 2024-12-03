"""
Classifies: CHEBI:26893 tetracyclic triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetracyclic_triterpenoid(smiles: str):
    """
    Determines if a molecule is a tetracyclic triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetracyclic triterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has a triterpenoid skeleton
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 4:
        return False, "Less than 4 rings found"

    # Check if there are at least 4 rings that are fused
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    fused_rings = 0
    for ring in atom_rings:
        for atom_idx in ring:
            if sum(atom_idx in other_ring for other_ring in atom_rings) > 1:
                fused_rings += 1
                break

    if fused_rings < 4:
        return False, "Less than 4 fused rings found"

    # Check if the molecule has a triterpenoid backbone
    # Triterpenoids generally have around 30 carbon atoms, but can have more
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if num_carbons < 27:  # Allowing for some flexibility
        return False, "Less than 27 carbon atoms found"

    return True, "Molecule is a tetracyclic triterpenoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26893',
                          'name': 'tetracyclic triterpenoid',
                          'definition': 'Any triterpenoid consisting of a '
                                        'tetracyclic skeleton.',
                          'parents': ['CHEBI:177333', 'CHEBI:36615']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 25,
    'num_false_positives': 19,
    'num_true_negatives': 1,
    'num_false_negatives': 0,
    'precision': 0.5681818181818182,
    'recall': 1.0,
    'f1': 0.7246376811594203,
    'accuracy': None}