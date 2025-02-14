"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: CHEBI:35980 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is formed by the condensation of tetradecanoic acid (myristic acid) with an alcohol or phenol.
    It contains a linear, unbranched acyl chain of 14 carbons (including the carbonyl carbon) connected via an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester functional group pattern
    # [C:1](=O)[O:2][C:3]
    ester_pattern = Chem.MolFromSmarts('[C:1](=O)[O:2][#6]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check if the acyl chain is tetradecanoyl
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon atom index
        ester_o_idx = match[1]     # Ester oxygen atom index
        acyl_chain_atoms = set()
        visited = set()

        # Traverse the acyl chain starting from the carbonyl carbon
        def traverse_acyl_chain(atom_idx, prev_idx=None):
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                return False  # Non-carbon atom in acyl chain
            if atom_idx in visited:
                return False  # Cycle detected
            visited.add(atom_idx)
            acyl_chain_atoms.add(atom_idx)
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() != prev_idx]
            
            # Remove the ester oxygen from neighbors
            neighbors = [idx for idx in neighbors if idx != ester_o_idx]

            # Check for branching
            carbon_neighbors = [idx for idx in neighbors if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
            if len(carbon_neighbors) > 1:
                return False  # Branching detected

            for neighbor_idx in carbon_neighbors:
                if neighbor_idx != prev_idx:
                    if not traverse_acyl_chain(neighbor_idx, atom_idx):
                        return False
            return True

        if traverse_acyl_chain(carbonyl_c_idx):
            # Check that the acyl chain length is 14 carbons
            if len(acyl_chain_atoms) == 14:
                return True, "Contains tetradecanoate ester group"
            else:
                continue  # Acyl chain length does not match
        else:
            continue  # Branching detected or invalid acyl chain

    return False, "No tetradecanoate ester groups found"