"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: CHEBI:18220 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible 1-benzopyran core pattern
    benzopyran_pattern = Chem.MolFromSmarts("O1C=C2C=CC=CC2=C1")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        # Try a more flexible pattern
        benzopyran_pattern = Chem.MolFromSmarts("O1C=C2C=CC=CC2=C1")
        if not mol.HasSubstructMatch(benzopyran_pattern):
            return False, "No 1-benzopyran core found"

    # Define the aryl substituent at position 3 pattern
    aryl_substituent_pattern = Chem.MolFromSmarts("O1C=C2C=CC=CC2=C1-c")
    if not mol.HasSubstructMatch(aryl_substituent_pattern):
        # Try a more flexible pattern
        aryl_substituent_pattern = Chem.MolFromSmarts("O1C=C2C=CC=CC2=C1-[c]")
        if not mol.HasSubstructMatch(aryl_substituent_pattern):
            return False, "No aryl substituent at position 3 found"

    # Check for the presence of at least one aromatic ring (aryl group)
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring)]
    if len(aromatic_rings) < 2:
        return False, "Not enough aromatic rings for isoflavonoid"

    # Check molecular weight - isoflavonoids typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for isoflavonoid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for isoflavonoid"
    if o_count < 1:
        return False, "Must have at least one oxygen (1-benzopyran core)"

    return True, "Contains 1-benzopyran core with aryl substituent at position 3"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18220',
                          'name': 'isoflavonoid',
                          'definition': 'Any 1-benzopyran with an aryl substituent at position 3. The term was originally restricted to natural products, but is now also used to describe semi-synthetic and fully synthetic compounds.',
                          'parents': ['CHEBI:50753', 'CHEBI:24835']},
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