"""
Classifies: CHEBI:65296 primary ammonium ion
"""
from rdkit import Chem

def is_primary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a primary ammonium ion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary ammonium ion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a primary ammonium ion group [NH3+]
    primary_ammonium_ion = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 1 and atom.GetTotalNumHs() == 3:
            primary_ammonium_ion = True
            break
    
    if not primary_ammonium_ion:
        return False, "No primary ammonium ion group [NH3+] found"
    
    # Check if the nitrogen is attached to a carbon (primary amino compound)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 1 and atom.GetTotalNumHs() == 3:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    return True, "Primary ammonium ion found"
    
    return False, "Nitrogen is not attached to a carbon (not a primary amino compound)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:65296',
                          'name': 'primary ammonium ion',
                          'definition': 'An ammonium ion derivative resulting '
                                        'from the protonation of the nitrogen '
                                        'atom of a primary amino compound. '
                                        'Major species at pH 7.3.',
                          'parents': ['CHEBI:25697', 'CHEBI:35274']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[02:25:29] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}