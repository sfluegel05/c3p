"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    Any fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated.
    This includes all fatty acyl-CoA(4-) molecules where the bond between the 11th and 12th carbons
    in the fatty acyl chain is either saturated (single bond) or absent (for chains shorter than 12 carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for fatty acyl-CoA thioester linkage
    # The fatty acyl chain is connected via a thioester linkage: C(=O)SCCNC(=O)
    # We will capture the carbonyl carbon of the fatty acyl chain
    fatty_acyl_pattern = Chem.MolFromSmarts('C(=O)SCCN')
    fatty_acyl_matches = mol.GetSubstructMatches(fatty_acyl_pattern)

    if not fatty_acyl_matches:
        return False, "No fatty acyl-CoA thioester linkage found"

    # Assume the first match is the fatty acyl chain
    fatty_acyl_match = fatty_acyl_matches[0]
    acyl_carbon_idx = fatty_acyl_match[0]  # Carbonyl carbon of the fatty acyl chain

    # Traverse the fatty acyl chain starting from the carbonyl carbon
    # We will build the acyl chain by moving away from the carbonyl carbon along carbons
    acyl_chain_atoms = []
    visited = set()

    def traverse_acyl_chain(atom_idx):
        """Recursively traverse the acyl chain starting from the given atom index."""
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom_idx in visited:
            return
        visited.add(atom_idx)

        if atom.GetAtomicNum() != 6:
            return  # Only consider carbon atoms

        acyl_chain_atoms.append(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
            if bond is None:
                continue
            if bond.IsInRing():
                continue  # Exclude ring structures
            if nbr_idx in visited:
                continue
            if nbr.GetAtomicNum() != 6:
                continue  # Only traverse to carbon atoms
            # Avoid going back towards the CoA moiety
            if nbr_idx == acyl_carbon_idx:
                continue
            traverse_acyl_chain(nbr_idx)

    # Start traversal from the carbon next to the carbonyl carbon
    carbonyl_carbon = mol.GetAtomWithIdx(acyl_carbon_idx)
    for neighbor in carbonyl_carbon.GetNeighbors():
        nbr_idx = neighbor.GetIdx()
        bond = mol.GetBondBetweenAtoms(acyl_carbon_idx, nbr_idx)
        if neighbor.GetAtomicNum() == 6 and bond.GetBondType() == rdchem.BondType.SINGLE:
            # Begin traversal
            traverse_acyl_chain(nbr_idx)
            break

    # Sort the acyl chain atoms based on their order from the carbonyl carbon
    # Here, acyl_chain_atoms[0] is the alpha carbon (adjacent to carbonyl carbon)
    # We need to ensure that the chain is correctly ordered
    # We can use a breadth-first search to get the correct order

    from collections import deque

    def get_chain_order(start_idx):
        """Get the ordered list of carbon indices in the acyl chain."""
        order = []
        queue = deque()
        visited = set()
        queue.append(start_idx)
        visited.add(start_idx)
        while queue:
            current_idx = queue.popleft()
            order.append(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)
            neighbors = [nbr for nbr in current_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            for nbr in neighbors:
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in visited and nbr_idx in acyl_chain_atoms:
                    queue.append(nbr_idx)
                    visited.add(nbr_idx)
        return order

    acyl_chain_ordered = get_chain_order(acyl_chain_atoms[0])

    # Count the number of carbons in the acyl chain
    num_carbons = len(acyl_chain_ordered)
    # For chains shorter than 12 carbons, we consider the C11-C12 bond as saturated (since it does not exist)
    if num_carbons < 11:
        return True, f"Fatty acyl chain has {num_carbons} carbons, C11-C12 bond is absent and considered saturated"

    # Get the atom indices for C11 and C12
    c11_idx = acyl_chain_ordered[10]  # 0-based indexing
    c12_idx = acyl_chain_ordered[11]

    # Get the bond between C11 and C12
    bond = mol.GetBondBetweenAtoms(c11_idx, c12_idx)
    if bond is None:
        return False, "No bond between C11 and C12"

    bond_type = bond.GetBondType()
    if bond_type == rdchem.BondType.SINGLE:
        return True, "C11-C12 bond is saturated (single bond)"
    elif bond_type == rdchem.BondType.DOUBLE or bond_type == rdchem.BondType.TRIPLE:
        return False, "C11-C12 bond is unsaturated (double or triple bond)"
    else:
        return False, f"C11-C12 bond is neither single nor double/triple (bond type: {bond_type})"