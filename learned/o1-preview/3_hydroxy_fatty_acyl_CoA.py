"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: CHEBI:134937 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is a coenzyme A derivative where the acyl group is a fatty acid chain 
    with a hydroxyl group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for Coenzyme A moiety
    coa_smarts = Chem.MolFromSmarts('NC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc23)[C@H](O)[C@@H]1OP(=O)(O)O')
    if not mol.HasSubstructMatch(coa_smarts):
        return False, "Coenzyme A moiety not found"

    # Define SMARTS pattern for thioester linkage to CoA
    thioester_smarts = Chem.MolFromSmarts('C(=O)SCCNC(=O)')
    if not mol.HasSubstructMatch(thioester_smarts):
        return False, "Thioester linkage to CoA not found"

    # Define SMARTS pattern for 3-hydroxy fatty acyl chain
    # This pattern looks for a chain of carbons connected via single bonds with a hydroxyl at the 3-position
    fatty_acyl_smarts = Chem.MolFromSmarts('C(=O)SC[C;H2][C;H](O)[C;H2]')
    matches = mol.GetSubstructMatches(fatty_acyl_smarts)
    if not matches:
        return False, "3-hydroxy fatty acyl chain not found"

    # Optionally, verify the length of the fatty acyl chain to ensure it's a fatty acid
    # Fatty acids typically have long hydrocarbon chains (usually at least 4 carbons)
    for match in matches:
        # The indices in the match correspond to the atoms in the SMARTS pattern
        # Match indices: [carbonyl C, S, first C, second C (with OH), third C]
        chain_atom_idx = match[-1]  # Last carbon in the pattern
        chain_length = 0
        visited_atoms = set(match)

        # Traverse the carbon chain beyond position 3
        atom = mol.GetAtomWithIdx(chain_atom_idx)
        while True:
            neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in visited_atoms]
            if not neighbors:
                break  # End of chain
            next_atom = neighbors[0]
            if next_atom.GetSymbol() != 'C':
                break  # Non-carbon atom, end traversal
            visited_atoms.add(next_atom.GetIdx())
            atom = next_atom
            chain_length += 1

        total_chain_length = chain_length + 3  # Including the first 3 carbons
        if total_chain_length >= 4:
            return True, "Contains 3-hydroxy fatty acyl-CoA structure"

    return False, "Does not match 3-hydroxy fatty acyl-CoA structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:134937',
        'name': '3-hydroxy fatty acyl-CoA',
        'definition': 'A hydroxy fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.',
        'parents': ['CHEBI:64479', 'CHEBI:83875']
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
    }
}