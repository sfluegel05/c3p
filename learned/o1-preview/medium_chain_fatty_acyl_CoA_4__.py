"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    A medium-chain fatty acyl-CoA(4-) is an acyl-CoA molecule with a fatty acyl chain
    of 6 to 12 carbons attached via a thioester bond to Coenzyme A, and with deprotonated
    phosphate groups resulting in a 4- charge at physiological pH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define Coenzyme A SMARTS pattern
    # This pattern represents the ADP and pantetheine portion of CoA
    # Using a simplified pattern to match the key features of CoA
    coa_smarts = Chem.MolFromSmarts("""
    [#8]-[#6]-[#1]-[#8]-[#3+](-[#8-])(-[#8-])-[#8]-[#6]-1-[#8]-[#6]-[#6]-[#8]-[#1]-1
    -[#7]-1-[#6]-[#7]-[#6]-2-[#7]-[#6]-[#7]-[#6]-2-[#7]-1
    """)

    if coa_smarts is None:
        return False, "Invalid CoA SMARTS pattern"

    if not mol.HasSubstructMatch(coa_smarts):
        return False, "Coenzyme A moiety not found"

    # Find the thioester bond (C(=O)-S-C)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S[C;!H0]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester bond not found"

    # Assume the first match is the acyl chain
    acyl_carbon_idx = thioester_matches[0][0]
    sulfur_idx = thioester_matches[0][2]

    # Traverse the acyl chain to count carbons
    visited = set()
    carbons_in_chain = 0
    stack = [mol.GetAtomWithIdx(acyl_carbon_idx)]

    while stack:
        atom = stack.pop()
        idx = atom.GetIdx()
        if idx in visited:
            continue
        visited.add(idx)
        if atom.GetAtomicNum() == 6:
            carbons_in_chain += 1
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Do not go back to sulfur atom
                if neighbor_idx == sulfur_idx:
                    continue
                stack.append(neighbor)

    # Subtract one to exclude the carbonyl carbon
    acyl_chain_length = carbons_in_chain - 1
    if acyl_chain_length < 6 or acyl_chain_length > 12:
        return False, f"Acyl chain length is {acyl_chain_length}, not in the range 6-12 for medium-chain fatty acids"

    # Check for deprotonated phosphate groups (4- charge)
    # Count the number of phosphate oxygens with negative charges
    negative_oxygen_count = sum(1 for atom in mol.GetAtoms()
                                if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if negative_oxygen_count < 4:
        return False, f"Found {negative_oxygen_count} negatively charged oxygens, expected at least 4 for deprotonated phosphates"

    return True, "Molecule is a medium-chain fatty acyl-CoA(4-) with appropriate acyl chain length and charge"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'medium-chain fatty acyl-CoA(4-)',
        'definition': 'An acyl-CoA oxoanion that results from deprotonation of the phosphate and diphosphate groups of any medium-chain fatty acyl-CoA; major species at pH 7.3.',
        'parents': []
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
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}