"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:57395 long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA results from the condensation of the thiol group of coenzyme A with
    the carboxy group of any long-chain (C13 to C22) fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the thioester linkage pattern S-C(=O)-C
    thioester_smarts = "[#16]-C(=O)-C"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pattern is None:
        return False, "Failed to construct thioester pattern"

    # Define adenine pattern (part of CoA)
    adenine_smarts = "n1cnc2c1ncnc2N"
    adenine_pattern = Chem.MolFromSmarts(adenine_smarts)
    if adenine_pattern is None:
        return False, "Failed to construct adenine pattern"

    # Check if adenine moiety is present (indicates CoA)
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Adenine moiety not found (CoA moiety missing)"

    # Find thioester linkage(s)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    for match in thioester_matches:
        sulfur_idx = match[0]
        carbonyl_c_idx = match[1]
        fatty_acyl_start_idx = match[2]

        # Traverse the acyl chain starting from fatty_acyl_start_idx
        visited = set()
        queue = [fatty_acyl_start_idx]
        carbon_count = 0

        while queue:
            current_idx = queue.pop(0)
            if current_idx in visited:
                continue
            visited.add(current_idx)

            atom = mol.GetAtomWithIdx(current_idx)
            if atom.GetAtomicNum() == 6:  # Carbon atom
                carbon_count += 1
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    # Exclude going back to the carbonyl carbon
                    if neighbor_idx not in visited and neighbor_idx != carbonyl_c_idx:
                        neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
                        if neighbor_atom.GetAtomicNum() == 6:
                            queue.append(neighbor_idx)
            else:
                # Stop traversal if atom is not carbon
                continue

        # Check if carbon count is within the long-chain fatty acid range (C13 to C22)
        if 13 <= carbon_count <= 22:
            return True, f"Contains long-chain fatty acyl group with {carbon_count} carbons"
        else:
            continue  # Check next thioester linkage if any

    return False, "No long-chain fatty acyl group of length 13-22 carbons found"