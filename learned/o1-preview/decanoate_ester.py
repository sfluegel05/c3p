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

    # Define ester group pattern with atom mapping
    ester_pattern = Chem.MolFromSmarts('[C:1](=O)[O:2][#6]')  # Ester group where :1 is carbonyl carbon, :2 is ester oxygen

    # Find ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check for decanoyl chain
    for match in ester_matches:
        if len(match) < 2:
            continue  # Skip if match does not have expected length

        carbonyl_c_idx = match[0]  # Carbonyl carbon atom index
        ester_o_idx = match[1]     # Ester oxygen atom index

        # Initialize variables for traversal
        visited = set()
        chain_length = 1  # Start counting from carbonyl carbon
        branching = False
        is_linear = True
        unsaturated = False

        # Traverse the acyl chain starting from carbonyl carbon
        def traverse_acyl_chain(atom_idx, prev_atom_idx):
            nonlocal chain_length, branching, is_linear, unsaturated
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)

            # Ensure atom is carbon
            if atom.GetAtomicNum() != 6:
                is_linear = False
                return

            # Ensure atom is sp3 hybridized (no double bonds)
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                unsaturated = True
                return

            # Check that the bond between prev_atom and current atom is single
            if prev_atom_idx != ester_o_idx:
                bond = mol.GetBondBetweenAtoms(prev_atom_idx, atom_idx)
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    unsaturated = True
                    return

            # Get neighbors excluding previous atom and ester oxygen atom
            neighbors = [nbr for nbr in atom.GetNeighbors() 
                         if nbr.GetIdx() != prev_atom_idx and nbr.GetIdx() != ester_o_idx]

            # Exclude non-carbon atoms
            neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 6]

            # If more than one neighbor, branching occurs
            if len(neighbors) > 1:
                branching = True
                is_linear = False
                return

            # Increment chain length
            if neighbors:
                chain_length += 1

            for neighbor in neighbors:
                if neighbor.GetIdx() not in visited:
                    traverse_acyl_chain(neighbor.GetIdx(), atom_idx)

        traverse_acyl_chain(carbonyl_c_idx, ester_o_idx)

        # Check if chain meets criteria
        if is_linear and not unsaturated and not branching and chain_length == 10:
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