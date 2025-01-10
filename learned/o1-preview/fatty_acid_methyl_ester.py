"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: CHEBI:XXXX fatty acid methyl ester
"""
from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is the carboxylic ester obtained by the formal condensation
    of a fatty acid with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the methyl ester SMARTS pattern
    methyl_ester_smarts = "[CX3](=O)[O][CH3]"
    methyl_ester_pattern = Chem.MolFromSmarts(methyl_ester_smarts)

    # Find all methyl ester groups in the molecule
    methyl_ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    num_methyl_esters = len(methyl_ester_matches)

    if num_methyl_esters == 0:
        return False, "No methyl ester groups found"

    # Check each methyl ester group to see if it's attached to a fatty acid chain
    for match in methyl_ester_matches:
        ester_carbon_idx = match[0]  # Carbonyl carbon atom index
        ester_oxygen_idx = match[1]  # Ester oxygen atom index

        # Get the carbonyl carbon atom
        carbonyl_carbon = mol.GetAtomWithIdx(ester_carbon_idx)

        # Find the atom attached to the carbonyl carbon that is not the ester oxygen
        neighbors = [nbr for nbr in carbonyl_carbon.GetNeighbors() if nbr.GetIdx() != ester_oxygen_idx]
        if len(neighbors) != 1:
            continue  # Move to the next match if not exactly one neighbor (should be the fatty acid chain)
        chain_start_atom = neighbors[0]

        # Traverse the chain to check if it's a fatty acid chain
        chain_atom_indices = set()
        atoms_to_visit = [chain_start_atom.GetIdx()]
        ring_found = False
        heteroatom_found = False

        while atoms_to_visit:
            current_idx = atoms_to_visit.pop()
            if current_idx in chain_atom_indices:
                continue
            chain_atom_indices.add(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)

            if current_atom.IsInRing():
                ring_found = True
                break  # Fatty acid chains should not contain rings

            atomic_num = current_atom.GetAtomicNum()
            if atomic_num != 6 and atomic_num != 1:  # Carbon or hydrogen
                heteroatom_found = True
                break  # Fatty acid chains should not contain heteroatoms

            # Add neighbor atoms to visit (excluding the ester carbon)
            for nbr in current_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx != carbonyl_carbon.GetIdx() and nbr_idx not in chain_atom_indices:
                    atoms_to_visit.append(nbr_idx)

        if ring_found:
            continue  # Not a fatty acid chain due to ring
        if heteroatom_found:
            continue  # Not a fatty acid chain due to heteroatom

        # Count the number of carbon atoms in the chain
        chain_length = sum(1 for idx in chain_atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if chain_length < 4:
            continue  # Chain too short to be a fatty acid

        # If we reach here, we have found a fatty acid methyl ester
        return True, "Molecule is a fatty acid methyl ester"

    # If no suitable fatty acid methyl ester groups are found
    return False, "No fatty acid methyl ester group found"