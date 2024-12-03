"""
Classifies: CHEBI:64674 1,2-diacyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_2_diacyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphoethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2-diacyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the phosphoethanolamine group
    phosphoethanolamine_pattern = Chem.MolFromSmarts("P(OC[C@H](O)COC(=O)*)OC(=O)CCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "Phosphoethanolamine group not found"

    # Check for the R-configuration at the glycerol backbone
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not any(center[1] == 'R' for center in chiral_centers if center[0] == 3):  # Assuming the chiral center is at index 3
        return False, "R-configuration at chiral center not found"

    return True, "1,2-diacyl-sn-glycero-3-phosphoethanolamine"

# Example usage
smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\C/C=C\C/C=C\CC)(OCCN)(O)=O"
result, reason = is_1_2_diacyl_sn_glycero_3_phosphoethanolamine(smiles)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64674',
                          'name': '1,2-diacyl-sn-glycero-3-phosphoethanolamine',
                          'definition': 'An optically active form of '
                                        'phosphatidylethanolamine having '
                                        'R-configuration.',
                          'parents': ['CHEBI:16038', 'CHEBI:36314']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'False Phosphoethanolamine group not found\n',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 62,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}