"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: unsaturated fatty acyl-CoA
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA is a coenzyme A linked via a thioester bond to an unsaturated fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define CoA substructure pattern (simplified)
    # CoA contains an adenosine diphosphate moiety and a pantetheine unit
    coa_smarts = "NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Error in CoA SMARTS pattern"

    # Check for CoA moiety
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"

    # Define thioester linkage pattern: C(=O)-S-
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if thioester_pattern is None:
        return False, "Error in thioester SMARTS pattern"

    # Find thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # For each thioester linkage, analyze the acyl chain
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon
        sulfur_idx = match[2] if len(match) >2 else match[1]  # Sulfur atom

        # Get the carbon attached to the carbonyl carbon (start of acyl chain)
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)

        # Identify the acyl chain starting atom (exclude carbonyl carbon)
        acyl_chain_atom = None
        for neighbor in carbonyl_c.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != sulfur_idx:
                acyl_chain_atom = neighbor
                break
        if acyl_chain_atom is None:
            continue  # No acyl chain found

        # Extract the acyl chain starting from acyl_chain_atom
        # Use BFS traversal to get the acyl chain atoms
        acyl_chain_atoms = set()
        to_visit = [acyl_chain_atom.GetIdx()]
        while to_visit:
            atom_idx = to_visit.pop()
            if atom_idx in acyl_chain_atoms:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            acyl_chain_atoms.add(atom_idx)
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                n_idx = neighbor.GetIdx()
                if n_idx in acyl_chain_atoms:
                    continue
                if neighbor.GetAtomicNum() != 6 and neighbor.GetAtomicNum() != 1:
                    continue  # Only consider carbon and hydrogen atoms
                if n_idx == carbonyl_c_idx:
                    continue  # Skip back to carbonyl carbon
                to_visit.append(n_idx)

        chain_length = len([idx for idx in acyl_chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() != 1]) + 1  # Include the carbonyl carbon

        if chain_length < 4:
            continue  # Chain too short to be a fatty acid

        # Check for presence of at least one C=C double bond in the acyl chain
        double_bond_found = False
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                if begin_idx in acyl_chain_atoms and end_idx in acyl_chain_atoms:
                    # Ensure the double bond is between carbons
                    if mol.GetAtomWithIdx(begin_idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(end_idx).GetAtomicNum() == 6:
                        double_bond_found = True
                        break

        if not double_bond_found:
            continue  # No double bond found in acyl chain

        return True, "Molecule is an unsaturated fatty acyl-CoA"

    return False, "No unsaturated fatty acid chain attached via thioester linkage found"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'unsaturated fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any unsaturated fatty acid.',
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
    'success': None,
    'best': None,
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