"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid is a fatty acid in which the hydrocarbon chain has one or more alkyl substituents.
    The fatty acyl chain is usually saturated and the substituent is often a methyl group; however, unsaturated BCFAs 
    and branches other than methyl are also possible.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[O;H1,-]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # For each carboxyl group, check the hydrocarbon chain
    for match in carboxylic_acid_matches:
        carboxyl_carbon_idx = match[0]
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

        # Start traversal from the carboxyl carbon to build the hydrocarbon chain
        visited = set()
        chain_atoms = []
        branches = 0
        stack = [(carboxyl_carbon_idx, None)]  # (current_atom_idx, previous_atom_idx)

        while stack:
            current_atom_idx, previous_atom_idx = stack.pop()
            if current_atom_idx in visited:
                continue
            visited.add(current_atom_idx)
            current_atom = mol.GetAtomWithIdx(current_atom_idx)

            # Include only carbon atoms in the chain
            if current_atom.GetAtomicNum() == 6:
                chain_atoms.append(current_atom_idx)
                neighbors = [nbr.GetIdx() for nbr in current_atom.GetNeighbors()]
                heavy_neighbors = [
                    nbr_idx for nbr_idx in neighbors
                    if mol.GetAtomWithIdx(nbr_idx).GetAtomicNum() > 1 and nbr_idx != previous_atom_idx
                ]

                # Check for branches (more than two heavy atom neighbors indicates branching)
                if len(heavy_neighbors) > 1:
                    branches += len(heavy_neighbors) - 1  # Exclude the main chain progression

                # Add neighbors to stack
                for nbr_idx in heavy_neighbors:
                    stack.append((nbr_idx, current_atom_idx))

        chain_length = len(chain_atoms)
        if chain_length < 4:
            continue  # Too short to be a fatty acid chain

        # Check for branches in the chain
        if branches > 0:
            return True, f"Contains carboxyl group with hydrocarbon chain of length {chain_length} and {branches} branch(es)"
        else:
            return False, "No branches detected in the hydrocarbon chain"

    # If no suitable chain with branches found
    return False, "No suitable hydrocarbon chain with branches found"