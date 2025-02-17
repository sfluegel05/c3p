"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: Fatty Aldehyde
Definition: An aldehyde formally arising from reduction of the carboxylic acid group of its corresponding fatty acid, having a carbonyl group at one end of the carbon chain.
"""

from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.

    A fatty aldehyde should have a long aliphatic carbon chain (we require at least 8 carbon atoms)
    and must contain at least one terminal aldehyde group. A terminal aldehyde group is defined here as:
      - An sp2 carbon that forms a C=O group,
      - Carries one hydrogen (as in an aldehyde),
      - Not part of a ring,
      - And attached to a carbon which is part of an aliphatic chain.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as fatty aldehyde, False otherwise.
        str: Explanation/reason for the classification.
    """
    # Parse the SMILES string to generate a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count total number of carbon atoms in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 8:
        return False, f"Not enough carbon atoms ({carbon_count} found; need at least 8) for a fatty chain"

    # Define a SMARTS pattern for an aldehyde group (the carbonyl carbon with one H)
    aldehyde_smarts = "[CX3H1](=O)"
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"

    # Check each aldehyde match for being terminal (i.e. at one end of a chain)
    terminal_aldehyde_found = False
    for match in aldehyde_matches:
        # match[0] will be the aldehyde carbon; match[1] is not explicitly provided by this SMARTS,
        # so we identify it through neighbors.
        aldehyde_c = mol.GetAtomWithIdx(match[0])
        # Ensure the aldehyde carbon is not in a ring (it should be terminal)
        if aldehyde_c.IsInRing():
            continue

        # Find the neighbor that is not oxygen. In an aldehyde, the carbonyl is double-bonded to O.
        neighbor_carbons = []
        for neighbor in aldehyde_c.GetNeighbors():
            # Oxygen neighbor is part of the C=O; skip it.
            if neighbor.GetAtomicNum() == 8:
                continue
            # We expect the other neighbor to be a carbon (attached to the fatty chain)
            if neighbor.GetAtomicNum() == 6:
                neighbor_carbons.append(neighbor)
        # For a terminal aldehyde group, there should be exactly one carbon neighbor.
        if len(neighbor_carbons) != 1:
            continue

        # Optionally, check that the attached carbon is part of a mostly linear aliphatic chain.
        # Here we simply check that it is an sp3 (or sp2) carbon and not in a ring.
        chain_carbon = neighbor_carbons[0]
        if chain_carbon.IsInRing():
            continue

        # If we reached here, we consider this aldehyde to be at the terminus of a fatty chain.
        terminal_aldehyde_found = True
        break

    if not terminal_aldehyde_found:
        return False, "No terminal aldehyde group found; aldehyde group(s) may be internal or in a ring"

    # Passed all checks: the molecule has a long carbon chain and a terminal aldehyde group.
    return True, "Molecule qualifies as a fatty aldehyde: has a long carbon chain with a terminal aldehyde group."


# Example test cases (uncomment to run)
# test_smiles = [
#     "[H]C(=CC=O)C(O)CCCCC",  # 4-hydroxynon-2-enal
#     "O=CCCCCCCCCC/C=C\\CCCCCCCC",  # 11Z-Eicosenal
#     "CCCCCCCCCCCC=O",  # dodecanal
#     "O=CCCCCCCCC/C=C\\CCCC",  # example dialdehyde possibility (adipaldehyde = O=CCCCCC=O if linear)
# ]
# for s in test_smiles:
#     result, reason = is_fatty_aldehyde(s)
#     print(f"SMILES: {s} -> {result}: {reason}")