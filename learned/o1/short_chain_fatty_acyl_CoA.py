"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:60940 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA results from the condensation of the thiol group of coenzyme A
    with the carboxy group of any short-chain fatty acid (2 to 5 carbons long).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for Coenzyme A moiety
    # Include key features: adenine ring, ribose, diphosphates, pantetheine arm
    # Use a simplified pattern to avoid over-specification
    coA_smarts_str = 'NC1=NC=NC2=C1N=CN2[C@H]3O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]4O[C@@H](N5C=NC=N5)[C@@H](O)[C@H]4OP(=O)(O)O)[C@H](O)[C@H]3OP(=O)(O)O'
    coA_pattern = Chem.MolFromSmarts(coA_smarts_str)

    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Coenzyme A moiety not found"

    # Define a SMARTS pattern for the thioester linkage
    thioester_smarts_str = 'C(=O)SC'
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts_str)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)

    if not thioester_matches:
        return False, "Thioester linkage not found"

    # For each thioester linkage, attempt to identify the acyl chain length
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Carbon of C=O
        sulfur_idx = match[2]      # Sulfur atom index

        # Get the acyl chain starting from the carbonyl carbon
        acyl_chain_atoms = []
        visited = set()
        stack = [carbonyl_c_idx]
        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            acyl_chain_atoms.append(atom)
            if atom_idx != carbonyl_c_idx:  # Do not include the carbonyl carbon again
                # Only traverse through carbons (exclude sulfur to prevent crossing into CoA moiety)
                if atom.GetAtomicNum() != 6:
                    continue
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Stop at the sulfur atom to prevent crossing into CoA
                if neighbor_idx == sulfur_idx:
                    continue
                if neighbor_idx not in visited:
                    stack.append(neighbor_idx)

        # Count the number of carbon atoms in the acyl chain, excluding the carbonyl carbon
        num_carbons = sum(1 for atom in acyl_chain_atoms if atom.GetAtomicNum() == 6) - 1  # Exclude carbonyl carbon

        if 1 <= num_carbons <= 4:
            return True, f"Found short-chain fatty acyl-CoA with acyl chain length {num_carbons +1} carbons"
        else:
            return False, f"Acyl chain length is {num_carbons +1}, not short-chain"

    return False, "Acyl chain not identified or not short-chain"

# Add metadata
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:60940',
        'name': 'short-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any short-chain fatty acid.',
        'parents': ['CHEBI:37554', 'CHEBI:57288']
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
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}