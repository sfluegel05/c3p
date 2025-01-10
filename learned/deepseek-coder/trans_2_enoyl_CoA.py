"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: CHEBI:28348 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    A trans-2-enoyl-CoA is characterized by a CoA moiety and a trans-2-enoyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts(
        "[*]SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    )
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Define the trans-2-enoyl group pattern (trans double bond at position 2, carbonyl at position 1)
    trans_2_enoyl_pattern = Chem.MolFromSmarts("[CX3]=[CX3H1]\\C=C\\C(=O)[SX2]")
    if not mol.HasSubstructMatch(trans_2_enoyl_pattern):
        return False, "No trans-2-enoyl group found"

    # Verify the thioester bond between the CoA thiol and the enoyl group
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester bond found between CoA and enoyl group"

    return True, "Contains CoA moiety and trans-2-enoyl group with a thioester bond"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:28348',
        'name': 'trans-2-enoyl-CoA',
        'definition': 'An unsaturated fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.',
        'parents': ['CHEBI:26348', 'CHEBI:76579']
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
    'accuracy': 0.9998521228585199
}