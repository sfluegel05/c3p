"""
Classifies: CHEBI:64612 1,2-diacyl-sn-glycero-3-phosphoethanolamine zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_1_2_diacyl_sn_glycero_3_phosphoethanolamine_zwitterion(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphoethanolamine zwitterion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2-diacyl-sn-glycero-3-phosphoethanolamine zwitterion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the zwitterion group [NH3+]
    zwitterion_group = Chem.MolFromSmarts("[NH3+]")
    if not mol.HasSubstructMatch(zwitterion_group):
        return False, "Zwitterion group [NH3+] not found"

    # Check for the presence of the phosphate group with [O-]
    phosphate_group = Chem.MolFromSmarts("P([O-])(=O)(OCC[NH3+])")
    if not mol.HasSubstructMatch(phosphate_group):
        return False, "Phosphate group with [O-] not found"

    # Check for the glycerol backbone with R-configuration
    glycerol_backbone = Chem.MolFromSmarts("[C@H](COC(=O)*)(COP([O-])(=O)OCC[NH3+])OC(=O)*")
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "Glycerol backbone with R-configuration not found"

    return True, "1,2-diacyl-sn-glycero-3-phosphoethanolamine zwitterion"

# Example usage:
# smiles = "[C@@H](COC(=O)*)(COP(OCC[NH3+])(=O)[O-])OC(=O)*"
# print(is_1_2_diacyl_sn_glycero_3_phosphoethanolamine_zwitterion(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64612',
                          'name': '1,2-diacyl-sn-glycero-3-phosphoethanolamine '
                                  'zwitterion',
                          'definition': 'An optically active '
                                        'phosphatidylethanolamine zwitterion '
                                        'having R-configuration',
                          'parents': ['CHEBI:57613']},
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
    'num_true_positives': 18,
    'num_false_positives': 0,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}