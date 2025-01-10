"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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

    # Add hydrogens to accurately count atoms
    mol = Chem.AddHs(mol)
    
    # Define Coenzyme A SMARTS pattern (simplified for matching)
    # This pattern represents the ADP and pantetheine portion of CoA
    coa_smarts = """
    O=P(O)([O-])OP(O)([O-])OC[C@H]1O[C@H](COP(O)([O-])=O)[C@@H](O)[C@H]1O
    n1cnc2c(ncnc12)N
    """
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Invalid CoA SMARTS pattern"

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"

    # Find the thioester bond (C(=O)-S-C)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester bond not found"

    # Locate the acyl chain attached to the thioester
    # Find the carbonyl carbon in the thioester bond
    thioester_matches = mol.GetSubstructMatch(thioester_pattern)
    acyl_carbon_idx = thioester_matches[0]
    acyl_carbon = mol.GetAtomWithIdx(acyl_carbon_idx)
    
    # Traverse the acyl chain to count carbons
    visited = set()
    carbons_in_chain = 0
    stack = [acyl_carbon]
    while stack:
        atom = stack.pop()
        idx = atom.GetIdx()
        if idx in visited:
            continue
        visited.add(idx)
        if atom.GetAtomicNum() == 6:
            carbons_in_chain += 1
            # Add neighboring atoms except those connected via double bonds to oxygen (to avoid counting carbonyl carbons multiple times)
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    continue
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    stack.append(neighbor)
                    
    # Subtract one to exclude the carbonyl carbon
    acyl_chain_length = carbons_in_chain - 1
    if acyl_chain_length < 6 or acyl_chain_length > 12:
        return False, f"Acyl chain length is {acyl_chain_length}, not in the range 6-12 for medium-chain fatty acids"

    # Check for deprotonated phosphate groups (4- charge)
    # Count the number of phosphate groups with negative charges
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])[O-]")
    phosphates = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphates) < 2:
        return False, "Phosphate groups not fully deprotonated to 4- charge state"

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
    'attempt': 0,
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