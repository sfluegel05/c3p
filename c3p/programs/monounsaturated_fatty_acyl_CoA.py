"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the CoA moiety
    coa_smarts = '[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12'
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"
    
    # Check for the presence of a fatty acyl chain with one double bond in the fatty acyl part
    double_bond_count = 0
    fatty_acyl_smarts = 'C(=O)SCCNC(=O)CCNC(=O)'
    fatty_acyl_pattern = Chem.MolFromSmarts(fatty_acyl_smarts)
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "Thioester linkage not found"
    
    fatty_acyl_part = mol.GetSubstructMatches(fatty_acyl_pattern)[0]
    fatty_acyl_atoms = set(fatty_acyl_part)
    
    # Find the fatty acyl chain which is the longest continuous chain starting from the carbonyl carbon
    carbonyl_carbon = None
    for atom_idx in fatty_acyl_part:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == 'C' and any(neighbor.GetSymbol() == 'O' for neighbor in atom.GetNeighbors()):
            carbonyl_carbon = atom
            break

    if carbonyl_carbon is None:
        return False, "Carbonyl carbon not found in fatty acyl chain"
    
    # Traverse the chain to find double bonds
    visited = set()
    stack = [carbonyl_carbon.GetIdx()]
    while stack:
        current_idx = stack.pop()
        if current_idx in visited:
            continue
        visited.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited and neighbor_idx not in fatty_acyl_atoms:
                stack.append(neighbor_idx)
                bond = mol.GetBondBetweenAtoms(current_idx, neighbor_idx)
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    double_bond_count += 1

    if double_bond_count != 1:
        return False, f"Expected 1 double bond in the fatty acyl chain, found {double_bond_count}"

    return True, "Monounsaturated fatty acyl-CoA"

# Example usage:
smiles = 'CCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12'
print(is_monounsaturated_fatty_acyl_CoA(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:139575',
                          'name': 'monounsaturated fatty acyl-CoA',
                          'definition': 'Any unsaturated fatty acyl-CoA in '
                                        'which the fatty acyl chain contains '
                                        'one carbon-carbon double bond.',
                          'parents': ['CHEBI:51006']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Monounsaturated fatty acyl-CoA')\n",
    'num_true_positives': 8,
    'num_false_positives': 2,
    'num_true_negatives': 8,
    'num_false_negatives': 2,
    'precision': 0.8,
    'recall': 0.8,
    'f1': 0.8000000000000002,
    'accuracy': None}