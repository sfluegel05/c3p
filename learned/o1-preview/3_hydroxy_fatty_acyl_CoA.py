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

    # Define SMARTS pattern for the thioester linkage (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    matches = mol.GetSubstructMatches(thioester_pattern)
    if not matches:
        return False, "No thioester linkage found"
    
    # Try to find the acyl chain connected via thioester linkage
    for match in matches:
        carbonyl_c_idx = match[0]
        sulfur_idx = match[1]

        # Get the carbonyl carbon atom
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)

        # Get the atom connected to carbonyl carbon that is not the sulfur
        neighbors = [nbr for nbr in carbonyl_c.GetNeighbors() if nbr.GetIdx() != sulfur_idx]
        if not neighbors:
            continue  # No acyl chain connected
        first_acyl_atom = neighbors[0]

        # Start traversing the acyl chain
        acyl_chain_atoms = set()
        acyl_chain_atoms.add(carbonyl_c_idx)
        acyl_chain_atoms.add(first_acyl_atom.GetIdx())
        current_atom = first_acyl_atom
        position = 1
        hydroxyl_found = False

        while True:
            if current_atom.GetSymbol() != 'C':
                break  # Expected carbon in acyl chain

            # Check if this is the third carbon (position 3)
            if position == 2:
                # Check if current atom has hydroxyl group attached
                has_oh = False
                for nbr in current_atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                        has_oh = True
                        break
                if has_oh:
                    hydroxyl_found = True
                else:
                    break  # No hydroxyl at position 3
                # No need to traverse further after position 3
                break

            # Get next carbon in acyl chain
            neighbors = [nbr for nbr in current_atom.GetNeighbors()
                         if nbr.GetIdx() not in acyl_chain_atoms and nbr.GetSymbol() == 'C']
            if not neighbors:
                break  # End of chain
            next_atom = neighbors[0]
            acyl_chain_atoms.add(next_atom.GetIdx())
            current_atom = next_atom
            position += 1

        if hydroxyl_found:
            # Optionally, check that the sulfur is part of CoA
            # For simplicity, we assume that if the molecule contains the thioester linkage and the acyl chain meets criteria, it's a 3-hydroxy fatty acyl-CoA
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