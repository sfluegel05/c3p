"""
Classifies: CHEBI:76579 triradylglycerol
"""
from rdkit import Chem

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for ester linkages at the three positions
    ester_patterns = [
        Chem.MolFromSmarts("OCC(O)COC(=O)"),
        Chem.MolFromSmarts("OCC(OC(=O))CO"),
        Chem.MolFromSmarts("OC(COC(=O))CO")
    ]
    
    ester_count = 0
    for pattern in ester_patterns:
        if mol.HasSubstructMatch(pattern):
            ester_count += 1

    if ester_count != 3:
        return False, "Not all three positions are esterified"

    return True, "Molecule is a triradylglycerol"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76579',
                          'name': 'triradylglycerol',
                          'definition': 'A glycerol compound having one of '
                                        'three possible substituent groups - '
                                        'either acyl, alkyl, or alk-1-enyl - '
                                        'at each of the three possible '
                                        'positions sn-1, sn-2 or sn-3. has '
                                        'functional parent glycerol '
                                        '(CHEBI:17754), children: triglyceride '
                                        '(CHEBI:17855). Parent: is_a '
                                        'glycerolipid (CHEBI:35741)',
                          'parents': ['CHEBI:35741']},
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
    'num_true_positives': 175,
    'num_false_positives': 14,
    'num_true_negatives': 6,
    'num_false_negatives': 0,
    'precision': 0.9259259259259259,
    'recall': 1.0,
    'f1': 0.9615384615384615,
    'accuracy': None}