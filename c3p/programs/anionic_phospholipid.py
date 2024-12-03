"""
Classifies: CHEBI:62643 anionic phospholipid
"""
from rdkit import Chem

def is_anionic_phospholipid(smiles: str):
    """
    Determines if a molecule is an anionic phospholipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anionic phospholipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a phosphate group with negative charge
    phosphate_group = Chem.MolFromSmarts('P([O-])(=O)O')
    if not mol.HasSubstructMatch(phosphate_group):
        return False, "No phosphate group with negative charge found"

    # Check for the presence of glycerol backbone
    glycerol_backbone = Chem.MolFromSmarts('OCC(O)CO')
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "No glycerol backbone found"

    # Check for the presence of long alkyl chains (indicative of lipids)
    long_alkyl_chain = Chem.MolFromSmarts('CCCCCCCCCCCCCCCC')
    if not mol.HasSubstructMatch(long_alkyl_chain):
        return False, "No long alkyl chains found"

    return True, "Molecule is an anionic phospholipid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:62643',
                          'name': 'anionic phospholipid',
                          'definition': 'Any organophosphate oxoanion that is '
                                        'a negatively charged phospholipid, '
                                        'e.g. phosphatidylserine(1-), '
                                        'phosphatidate(2-), '
                                        'phosphatidylglycerol(1-).',
                          'parents': ['CHEBI:58945']},
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
    'num_true_positives': 21,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 78,
    'precision': 0.9130434782608695,
    'recall': 0.21212121212121213,
    'f1': 0.34426229508196726,
    'accuracy': None}