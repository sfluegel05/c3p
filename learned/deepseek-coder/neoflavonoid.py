"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: CHEBI:XXXXX neoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible 1-benzopyran (chromene) core structure
    benzopyran_pattern = Chem.MolFromSmarts("[O]1[C]2=[C][C]=[C][C]=[C]2[C]1")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 1-benzopyran (chromene) core found"

    # Define the aryl substituent pattern
    aryl_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1")

    # Check if the aryl group is attached to the 4th position of the benzopyran core
    matches = mol.GetSubstructMatches(benzopyran_pattern)
    for match in matches:
        # The 4th position in the benzopyran core is the 3rd atom in the match (0-based index)
        core_atom_idx = match[2]
        core_atom = mol.GetAtomWithIdx(core_atom_idx)
        for neighbor in core_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                neighbor_mol = Chem.RWMol(mol)
                neighbor_mol.RemoveAtom(neighbor.GetIdx())
                if neighbor_mol.HasSubstructMatch(aryl_pattern):
                    return True, "Contains 1-benzopyran core with aryl substituent at position 4"

    return False, "No aryl substituent found at position 4 of the 1-benzopyran core"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:XXXXX',
                          'name': 'neoflavonoid',
                          'definition': 'Any 1-benzopyran with an aryl substituent at position 4.',
                          'parents': ['CHEBI:XXXXX', 'CHEBI:XXXXX']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}