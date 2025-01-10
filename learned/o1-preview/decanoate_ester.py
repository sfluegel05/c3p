"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:36027 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester is a fatty acid ester resulting from the formal condensation of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester group pattern: carbonyl carbon single bonded to oxygen
    ester_pattern = Chem.MolFromSmarts('[$(C(=O)[O;!$([O-])])]')  # Matches ester functional group

    # Find ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check for decanoyl chain
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon atom index
        ester_o_idx = match[1]     # Ester oxygen atom index

        # Initialize variables
        visited = set()
        chain_length = 1  # Start counting from carbonyl carbon
        branching = False
        is_linear = True

        # Traverse the acyl chain starting from carbonyl carbon
        def traverse_acyl_chain(atom_idx):
            nonlocal chain_length, branching, is_linear
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)

            # Only consider carbon atoms
            if atom.GetAtomicNum() != 6:
                return

            neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() != ester_o_idx and nbr.GetIdx() not in visited]

            # Check for branching
            num_carbon_neighbors = sum(1 for nbr in neighbors if nbr.GetAtomicNum() == 6)
            if num_carbon_neighbors > 1:
                branching = True
                is_linear = False
                return  # Stop traversal if branching occurs

            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 6:
                    chain_length += 1
                    traverse_acyl_chain(neighbor.GetIdx())

        traverse_acyl_chain(carbonyl_c_idx)

        # Check if chain is linear and has exactly 10 carbons
        if not branching and chain_length == 10:
            return True, "Contains decanoate ester group"

    return False, "No decanoate ester group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36027',
        'name': 'decanoate ester',
        'definition': 'A fatty acid ester resulting from the formal condensation of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.',
        'parents': ['CHEBI:35620']
    }
}