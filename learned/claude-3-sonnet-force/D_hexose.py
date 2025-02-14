"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:16646 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose that has D-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains exactly 6 carbon atoms and 6 oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count != 6 or o_count != 6:
        return False, "Not a hexose (does not contain 6 carbons and 6 oxygens)"

    # Check for D-configuration at position 5 (carbon with 3 substituents)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnspec=True)
    has_d_config_at_5 = False
    for chiral_center in chiral_centers:
        atom = mol.GetAtomWithIdx(chiral_center[0])
        if atom.GetDegree() == 4 and atom.GetTotalNumHs() == 0:  # Carbon with 3 substituents
            chiral_tag = Chem.FindMolChiralCenterStereoTag(mol, chiral_center[0], chiral_center[1])
            if chiral_tag == Chem.MDL_R_CHIRALITY:
                has_d_config_at_5 = True
                break

    if not has_d_config_at_5:
        return False, "Does not have D-configuration at position 5"

    return True, "Contains 6 carbons, 6 oxygens, and has D-configuration at position 5"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:16646',
        'name': 'D-hexose',
        'definition': 'A hexose that has D-configuration at position 5.',
        'parents': ['CHEBI:18237', 'CHEBI:24051']
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
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 184,
    'num_false_positives': 1,
    'num_true_negatives': 182426,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.9947916666666666,
    'recall': 0.994840694978992,
    'f1': 0.9948161631216901,
    'accuracy': 0.9998505063367151
}