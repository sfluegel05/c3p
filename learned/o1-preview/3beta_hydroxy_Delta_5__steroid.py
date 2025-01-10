"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid

Definition: 'Any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6.'
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid is any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6 (Delta5).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add hydrogens to the molecule to ensure proper matching of stereochemistry
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates to assign stereochemistry if not present
    AllChem.EmbedMolecule(mol, maxAttempts=50, enforceChirality=True)
    Chem.AssignAtomChiralTagsFromStructure(mol)
    
    # Define steroid nucleus pattern (cyclopentanophenanthrene core with correct ring fusion)
    steroid_nucleus = Chem.MolFromSmarts("""
    [C@H]1CC[C@H]2[C@@H]3CC=C4[C@@H]([C@]4(CC[C@@H]3CC2)C)[C@@H]1C
    """)
    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "No steroid nucleus with correct stereochemistry found"
    
    # Define 3beta-hydroxy group pattern with beta orientation
    hydroxy_3beta = Chem.MolFromSmarts("""
    [C@@H]([C@@H]1CC[C@H]2[C@@H]3CC=C4[C@@H]([C@]4(CC[C@@H]3CC2)C)[C@@H]1C)O
    """)
    if not mol.HasSubstructMatch(hydroxy_3beta):
        return False, "No 3beta-hydroxy group with beta orientation found"
    
    # Define Delta(5) double bond between positions 5 and 6
    delta5_double_bond = Chem.MolFromSmarts("""
    [C@H]1CC=C[C@@H](C)CC1
    """)
    if not mol.HasSubstructMatch(delta5_double_bond):
        return False, "No Delta(5) double bond between positions 5 and 6 found"
    
    return True, "Molecule is a 3beta-hydroxy-Delta(5)-steroid"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': '3beta-hydroxy-Delta(5)-steroid',
        'definition': 'Any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6.'
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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}