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

    # More specific trans-2-enoyl pattern:
    # - Must have exactly one trans double bond at position 2
    # - Must have exactly one carbonyl at position 1
    # - Must be connected to CoA via thioester bond
    trans_2_enoyl_pattern = Chem.MolFromSmarts("[CX3H1,CX4H2]\\C=C\\C(=O)[SX2]")
    matches = mol.GetSubstructMatches(trans_2_enoyl_pattern)
    
    if not matches:
        return False, "No trans-2-enoyl group found"
    
    # Verify there's exactly one trans-2-enoyl group
    if len(matches) != 1:
        return False, f"Found {len(matches)} potential enoyl groups, need exactly 1"

    # Verify the double bond is truly trans
    bond = mol.GetBondBetweenAtoms(matches[0][1], matches[0][2])
    if not bond or not bond.GetStereo() == Chem.rdchem.BondStereo.STEREOE:
        return False, "Double bond is not in trans configuration"

    # Verify the thioester bond is connected to CoA
    thioester_atom = matches[0][-1]  # The sulfur atom
    for neighbor in mol.GetAtomWithIdx(thioester_atom).GetNeighbors():
        if neighbor.GetSymbol() == "C" and neighbor.GetDegree() == 3:  # CoA connection
            return True, "Contains CoA moiety and trans-2-enoyl group with a thioester bond"

    return False, "Thioester bond not properly connected to CoA moiety"


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