"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: CHEBI:134937 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # Define SMARTS pattern for Coenzyme A thioester linkage
    # The pattern includes the CoA moiety attached via a thioester linkage to an acyl chain
    coa_thioester_pattern = Chem.MolFromSmarts("S[C](=O)[C]")  # Simplified pattern for thioester

    if not mol.HasSubstructMatch(coa_thioester_pattern):
        return False, "No CoA thioester linkage found"

    # Define SMARTS pattern for 3-hydroxy acyl chain
    # This pattern looks for S-C(=O)-C-C-C(O)- chain
    hydroxy_acyl_pattern = Chem.MolFromSmarts("S[C](=O)[C][C][C](O)")  # 3rd carbon has OH

    if not mol.HasSubstructMatch(hydroxy_acyl_pattern):
        return False, "No 3-hydroxy fatty acyl chain found"

    # Optionally, we can check that the acyl chain is of sufficient length (e.g., more than 4 carbons)
    # Count the number of carbons in the acyl chain
    # Starting from the carbonyl carbon, traverse the chain

    # Find the matches for the thioester linkage
    matches = mol.GetSubstructMatches(coa_thioester_pattern)
    found = False
    for match in matches:
        sulfur_idx = match[0]
        carbonyl_c_idx = match[1]
        first_c_idx = match[2]

        # Now, traverse the acyl chain starting from first_c_idx
        acyl_chain_atoms = [carbonyl_c_idx, first_c_idx]
        current_atom_idx = first_c_idx
        chain_length = 1  # Start counting from the first carbon after carbonyl

        # We will build a list of atoms in the acyl chain
        while True:
            atom = mol.GetAtomWithIdx(current_atom_idx)
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() not in acyl_chain_atoms and nbr.GetAtomicNum() == 6]  # Only carbon atoms
            if len(neighbors) == 0:
                break  # End of chain
            elif len(neighbors) > 1:
                # Branching occurs, but fatty acyl chains are unbranched
                return False, "Branching found in acyl chain"
            else:
                next_atom_idx = neighbors[0]
                acyl_chain_atoms.append(next_atom_idx)
                current_atom_idx = next_atom_idx
                chain_length += 1

                # Check if current atom is the third carbon (C3)
                if chain_length == 3:
                    # Check if this carbon has a hydroxyl group
                    c3_atom = mol.GetAtomWithIdx(current_atom_idx)
                    has_hydroxy = False
                    for nbr in c3_atom.GetNeighbors():
                        if nbr.GetAtomicNum() == 8:  # Oxygen
                            if nbr.GetDegree() == 1:  # Only connected to one atom
                                has_hydroxy = True
                                break
                    if not has_hydroxy:
                        return False, "No hydroxyl group at position 3 of acyl chain"
                    else:
                        found = True
                        break
        if found:
            return True, "Contains 3-hydroxy fatty acyl-CoA structure"

    return False, "No 3-hydroxy fatty acyl-CoA structure found"

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
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 100,
    'num_false_positives': 5,
    'num_true_negatives': 182500,
    'num_false_negatives': 10,
    'num_negatives': None,
    'precision': 0.9523809523809523,
    'recall': 0.9090909090909091,
    'f1': 0.9302325581395349,
    'accuracy': 0.9999172758656234
}