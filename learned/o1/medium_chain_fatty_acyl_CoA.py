"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:60903 medium-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA results from the condensation of coenzyme A with a medium-chain fatty acid.
    Medium-chain fatty acids are aliphatic chains with 6 to 13 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thioester linkage pattern: C(=O)S
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # For each thioester linkage, attempt to find the acyl chain
    for match in thioester_matches:
        carbonyl_c_idx = match[0]
        sulfur_idx = match[1]

        # Initialize variables for traversal
        visited = set()
        carbons_in_chain = set()
        queue = []

        # Start from the carbonyl carbon, go along the acyl chain away from the sulfur
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        for neighbor in carbonyl_c.GetNeighbors():
            if neighbor.GetIdx() != sulfur_idx:
                queue.append(neighbor)

        contains_aromatic = False

        while queue:
            atom = queue.pop()
            atom_idx = atom.GetIdx()

            if atom_idx in visited:
                continue
            visited.add(atom_idx)

            if atom.GetAtomicNum() == 6:
                carbons_in_chain.add(atom_idx)
                if atom.GetIsAromatic():
                    contains_aromatic = True
            else:
                # If we encounter non-carbon atoms (other than O in functional groups), we may be entering CoA or other moieties
                pass

            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in visited:
                    continue
                # Avoid going back to the carbonyl carbon or towards the sulfur
                if neighbor_idx != carbonyl_c_idx and neighbor_idx != sulfur_idx:
                    bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
                    queue.append(neighbor)

        # Exclude aromatic acyl chains
        if contains_aromatic:
            continue  # Try next thioester linkage

        # Count the number of carbons in the acyl chain
        carbon_count = len(carbons_in_chain)

        # Medium-chain fatty acids have 6 to 13 carbons
        if 6 <= carbon_count <= 13:
            return True, f"Contains fatty acyl chain with {carbon_count} carbons (medium-chain)"
        else:
            return False, f"Fatty acyl chain with {carbon_count} carbons is not medium-chain"

    # If no suitable acyl chain found
    return False, "No suitable fatty acyl chain found"