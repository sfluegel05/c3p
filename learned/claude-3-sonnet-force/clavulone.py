"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: CHEBI:51721 clavulone
A class of esterified prostanoids obtained from marine corals.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    Clavulones are a class of esterified prostanoids obtained from marine corals,
    typically containing a cyclopentenone core, trans double bonds, halogen substituents,
    ester groups, and long carbon chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyclopentenone core
    cyclopentenone_pattern = Chem.MolFromSmarts("[C@H]1=C([C@@]2(C=C1)O[C@@H]2[H])O")
    if not mol.HasSubstructMatch(cyclopentenone_pattern):
        return False, "No cyclopentenone core found"
    
    # Look for trans double bonds in side chains
    trans_double_bond_pattern = Chem.MolFromSmarts("/C=C/CCCC")
    if not mol.HasSubstructMatch(trans_double_bond_pattern):
        return False, "No trans double bonds in side chains"
    
    # Look for halogen substituents (Cl, Br, I)
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I]")
    if not mol.HasSubstructMatch(halogen_pattern):
        return False, "No halogen substituent found"
    
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"
    
    # Look for long carbon chains (at least 6 carbons)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain found"
    
    return True, "Contains cyclopentenone core, trans double bonds, halogen substituent, ester groups, and long carbon chains"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:51721',
        'name': 'clavulone',
        'definition': 'A class of esterified prostanoids obtained from marine corals.',
        'parents': ['CHEBI:23995', 'CHEBI:39679']
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
    'message': 'Your program successfully classified the provided examples and achieved an F1 score above the threshold. Well done!',
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 124,
    'num_false_positives': 8,
    'num_true_negatives': 182405,
    'num_false_negatives': 46,
    'num_negatives': None,
    'precision': 0.9394,
    'recall': 0.7295,
    'f1': 0.8226,
    'accuracy': 0.9998
}