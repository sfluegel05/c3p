"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a fatty acid containing only methyl branches.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group and identify the carboxyl carbon
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxylic_matches:
        return False, "No carboxylic acid group found"

    # Assume the first carboxylic acid group if multiple are present
    carboxyl_c_idx = carboxylic_matches[0][0]
    carboxyl_c_atom = mol.GetAtomWithIdx(carboxyl_c_idx)

    # Function to get the main carbon chain starting from the carboxyl carbon
    def get_main_chain(atom_idx, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom_idx in visited:
            return []
        visited.add(atom_idx)

        if atom.GetAtomicNum() != 6:
            # Stop at heteroatoms (e.g., oxygen in carboxylic acid)
            return []

        chain = [atom_idx]
        for neighbor in atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 6 and nbr_idx not in visited:
                chain.extend(get_main_chain(nbr_idx, visited))
        return chain

    # Get the main chain
    visited_atoms = set()
    main_chain = get_main_chain(carboxyl_c_idx, visited_atoms)

    if len(main_chain) < 2:
        return False, "Main carbon chain is too short"

    # Check for branches off the main chain
    methyl_branch_count = 0
    for atom_idx in main_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx not in main_chain:
                # Found a branch
                branch_atom = mol.GetAtomWithIdx(nbr_idx)
                if branch_atom.GetAtomicNum() != 6:
                    return False, "Branch contains heteroatoms"
                # Check if the branch is a methyl group (terminal carbon)
                # Should have only one heavy atom neighbor (the main chain carbon)
                if len([nbr for nbr in branch_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]) > 1:
                    return False, "Branch larger than methyl group found"
                else:
                    methyl_branch_count += 1

    # Ensure there is at least one methyl branch
    if methyl_branch_count == 0:
        return False, "No methyl branches found"

    return True, "Molecule is a methyl-branched fatty acid with only methyl branches"

__metadata__ = {
    'chemical_class': {
        'name': 'methyl-branched fatty acid',
        'definition': 'Any branched-chain fatty acid containing methyl branches only.'
    },
    'message': None,
    'success': True,
    'error': '',
}