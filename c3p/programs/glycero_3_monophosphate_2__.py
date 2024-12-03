"""
Classifies: CHEBI:78185 glycero-3-monophosphate(2-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_glycero_3_monophosphate_2__(smiles: str):
    """
    Determines if a molecule is a glycero-3-monophosphate(2-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycero-3-monophosphate(2-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the phosphate group at sn-3 position
    phosphate_group = Chem.MolFromSmarts('COP([O-])([O-])=O')
    if not mol.HasSubstructMatch(phosphate_group):
        return False, "No phosphate group at sn-3 position found"

    # Check for the glycerol backbone
    glycerol_backbone = Chem.MolFromSmarts('O[C@@H](CO[*])[*]')
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "No glycerol backbone found"

    # Check for the presence of acyl, alkyl, or alkenyl groups
    acyl_group = Chem.MolFromSmarts('C(=O)O[*]')
    alkyl_group = Chem.MolFromSmarts('C[*]')
    alkenyl_group = Chem.MolFromSmarts('C=C[*]')

    has_acyl_group = mol.HasSubstructMatch(acyl_group)
    has_alkyl_group = mol.HasSubstructMatch(alkyl_group)
    has_alkenyl_group = mol.HasSubstructMatch(alkenyl_group)

    if not (has_acyl_group or has_alkyl_group or has_alkenyl_group):
        return False, "No acyl, alkyl, or alkenyl groups found"

    # Check for ester, ether, or vinyl linkages
    ester_linkage = Chem.MolFromSmarts('C(=O)O[*]')
    ether_linkage = Chem.MolFromSmarts('CO[*]')
    vinyl_linkage = Chem.MolFromSmarts('C=C[*]')

    has_ester_linkage = mol.HasSubstructMatch(ester_linkage)
    has_ether_linkage = mol.HasSubstructMatch(ether_linkage)
    has_vinyl_linkage = mol.HasSubstructMatch(vinyl_linkage)

    if not (has_ester_linkage or has_ether_linkage or has_vinyl_linkage):
        return False, "No ester, ether, or vinyl linkages found"

    return True, "Molecule is a glycero-3-monophosphate(2-)"

# Example usage:
smiles = 'CCCCCCCC\C=C/CCCCCCCC(=O)O[C@H](COC(=O)CCCCCCC\C=C/C\C=C/C\C=C/CC)COP([O-])([O-])=O'
print(is_glycero_3_monophosphate_2__(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:78185',
                          'name': 'glycero-3-monophosphate(2-)',
                          'definition': 'An anionic phospholipid having a '
                                        'phosphate group at sn-1 or sn-3 '
                                        'position of the glycerol backbone, '
                                        'and with a combination of one or two '
                                        'acyl groups, alkyl groups, or alkenyl '
                                        'groups  attached at the sn-1, sn-2, '
                                        'or sn-3 positions through ester, '
                                        'ether or vinyl linkages respectively.',
                          'parents': ['CHEBI:62643']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Molecule is a glycero-3-monophosphate(2-)')\n",
    'num_true_positives': 11,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}