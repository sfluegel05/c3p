"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:84947 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    An 11,12-saturated fatty acyl-CoA(4-) is a fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA head group
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12)C(=O)NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA head group"
    
    # Look for fatty acyl chain (at least 12 carbons)
    fatty_acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])CC[CX4,CX3]([H])([H])[CX4,CX3]([H])([H])[CX4,CX3]([H])([H])[CX4,CX3]([H])([H])[CX4,CX3]([H])([H])[CX4,CX3]([H])([H])[CX4,CX3]([H])([H])[CX4,CX3]([H])([H])[CX4,CX3]([H])([H])[CX4,CX3]([H])([H])")
    fatty_acyl_matches = mol.GetSubstructMatches(fatty_acyl_pattern)
    if not fatty_acyl_matches:
        return False, "Missing fatty acyl chain (at least 12 carbons)"
    
    # Check for 11-12 saturated bond
    for match in fatty_acyl_matches:
        chain = mol.GetAtomWithIdx(match[0]).GetNeighbors()
        for i in range(10, len(chain)):
            atom1 = chain[i]
            atom2 = chain[i+1]
            bond = mol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx())
            if bond.GetBondType() == Chem.BondType.SINGLE and i == 10:
                return True, "Contains an 11-12 saturated bond in the fatty acyl chain"
    
    return False, "No 11-12 saturated bond found in the fatty acyl chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:84947',
        'name': '11,12-saturated fatty acyl-CoA(4-)',
        'definition': 'Any fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated. This is needed for a reaction together with CHEBI:84947.',
        'parents': ['CHEBI:35616', 'CHEBI:35618']
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
    'stdout': None,
    'num_true_positives': 54,
    'num_false_positives': 0,
    'num_true_negatives': 182434,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 1.0
}