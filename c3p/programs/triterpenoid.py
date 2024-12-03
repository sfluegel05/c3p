"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has 30 carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
    if len(carbon_atoms) < 30:
        return False, "Less than 30 carbon atoms found"

    # Check if the molecule has ring structures
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No ring structures found"

    # Check for common triterpenoid modifications
    num_oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    if num_oxygen_atoms < 1:
        return False, "No oxygen atoms found, which are typical in triterpenoids"

    # Check if the molecule is derived from a triterpene
    # This part is heuristic and may need refinement
    # Assuming triterpenoids have at least one ring structure with more than 5 carbon atoms
    has_large_ring = any(len(ring) > 5 for ring in ring_info.AtomRings())
    if not has_large_ring:
        return False, "No large ring structures found"

    return True, "Molecule meets basic criteria for triterpenoids"

# Example usage:
# smiles = "O(C1C(C2C3(C4(C3)C(C5(C(CC4)(C(CC5)C(CCCC(C)C)C)C)C)CC2)CC1)(C)C)C(=O)/C=C/C6=CC(OC)=C(O)C=C6"
# print(is_triterpenoid(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36615',
                          'name': 'triterpenoid',
                          'definition': 'Any terpenoid derived from a '
                                        'triterpene. The term includes '
                                        'compounds in which the C30 skeleton '
                                        'of the parent triterpene has been '
                                        'rearranged or modified by the removal '
                                        'of one or more skeletal atoms '
                                        '(generally methyl groups).',
                          'parents': ['CHEBI:26873']},
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
    'num_true_positives': 211,
    'num_false_positives': 6,
    'num_true_negatives': 14,
    'num_false_negatives': 28,
    'precision': 0.9723502304147466,
    'recall': 0.8828451882845189,
    'f1': 0.9254385964912281,
    'accuracy': None}