"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule contains a glycerol backbone
    glycerol_smiles = "OCC(O)CO"
    glycerol_mol = Chem.MolFromSmiles(glycerol_smiles)
    if not mol.HasSubstructMatch(glycerol_mol):
        return False, "No glycerol backbone found"
    
    # Check for one acyl, alkyl or alk-1-enyl substituent
    acyl_pattern = Chem.MolFromSmarts("C(=O)[C;!R]")
    alkyl_pattern = Chem.MolFromSmarts("C[C;!R]")
    alk1enyl_pattern = Chem.MolFromSmarts("C=C[C;!R]")
    
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)
    alk1enyl_matches = mol.GetSubstructMatches(alk1enyl_pattern)
    
    total_substituents = len(acyl_matches) + len(alkyl_matches) + len(alk1enyl_matches)
    
    if total_substituents != 1:
        return False, f"Expected 1 acyl, alkyl, or alk-1-enyl substituent, found {total_substituents}"
    
    return True, "Molecule is a monoradylglycerol"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76575',
                          'name': 'monoradylglycerol',
                          'definition': 'Any lipid that is glycerol bearing a '
                                        'single acyl, alkyl or alk-1-enyl '
                                        'substituent at an unspecified '
                                        'position.',
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 16,
    'num_false_negatives': 16,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}