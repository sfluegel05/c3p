"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:57395 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA results from the condensation of the thiol group of coenzyme A with
    the carboxy group of any long-chain (C13 to C22) fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the Coenzyme A (CoA) pattern
    coa_smarts = """
    C[C@@H](O)[C@@H](COP(O)(=O)OCC1OC(O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O)n1cnc2c(N)ncnc12
    """
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Failed to construct CoA pattern"

    # Check for CoA substructure
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"

    # Define thioester linkage pattern (S-C(=O)-)
    thioester_pattern = Chem.MolFromSmarts("SC(=O)C")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Find the fatty acyl chain attached via thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Identify the fatty acyl chain
    for match in thioester_matches:
        sulfur_idx = match[0]
        # Traverse the chain from the carbonyl carbon
        carbonyl_c_idx = match[1]
        fatty_acyl_chain = Chem.FragmentOnBonds(mol, [carbonyl_c_idx], addDummies=True)
        # Extract the fatty acyl fragment
        fragments = Chem.GetMolFrags(fatty_acyl_chain, asMols=True, sanitizeFrags=False)
        for frag in fragments:
            atom_nums = [atom.GetAtomicNum() for atom in frag.GetAtoms()]
            if 6 in atom_nums and 1 in atom_nums:  # Contains carbon and hydrogen
                # Count number of carbons excluding the carbonyl carbon
                c_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
                c_count -= 1  # Exclude the carbonyl carbon
                if 13 <= c_count <= 22:
                    return True, f"Contains long-chain fatty acyl group with {c_count} carbons"
                else:
                    return False, f"Fatty acyl chain length is {c_count} carbons, not in range 13-22"

    return False, "Failed to identify long-chain fatty acyl chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:57395',
        'name': 'long-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.',
        'parents': ['CHEBI:37554', 'CHEBI:64479']
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