"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is an unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify thioester linkage: C(=O)-S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Assume first match is the thioester linkage
    carbonyl_c_idx = thioester_matches[0][0]
    sulfur_idx = thioester_matches[0][2]
    sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)

    # Check for CoA moiety connected to sulfur atom
    # From sulfur atom, check for at least two phosphate groups and an adenine ring

    # Function to perform BFS traversal from sulfur atom
    def traverse_from_sulfur(sulfur_idx):
        visited = set()
        queue = [sulfur_idx]
        phosphate_count = 0
        adenine_found = False

        while queue:
            atom_idx = queue.pop(0)
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)

            # Check for phosphate group [P](=O)(O)(O)
            if atom.GetAtomicNum() == 15:
                # Phosphorus atom found, check connections
                o_count = 0
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8:
                        o_count += 1
                if o_count >= 3:
                    phosphate_count += 1

            # Check for adenine ring
            # Adenine has a purine ring system with specific nitrogen patterns
            # SMARTS pattern for adenine
            adenine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2N")
            if mol.HasSubstructMatch(adenine_pattern):
                adenine_found = True

            # Add neighbors to queue
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    queue.append(neighbor_idx)

        return phosphate_count, adenine_found

    phosphate_count, adenine_found = traverse_from_sulfur(sulfur_idx)
    if phosphate_count < 2:
        return False, f"Less than 2 phosphate groups found ({phosphate_count})"
    if not adenine_found:
        return False, "Adenine ring of CoA not found"

    # Now check for double bond between positions 2 and 3 in the acyl chain
    # Starting from carbonyl carbon, traverse the acyl chain

    # Get carbonyl carbon atom
    carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
    # Find the alpha carbon (next carbon in acyl chain)
    alpha_c = None
    for neighbor in carbonyl_c.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != sulfur_idx:
            alpha_c = neighbor
            break
    if alpha_c is None:
        return False, "Alpha carbon in acyl chain not found"

    # Find beta carbon (next carbon after alpha)
    beta_c = None
    for neighbor in alpha_c.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != carbonyl_c_idx:
            beta_c = neighbor
            break
    if beta_c is None:
        return False, "Beta carbon in acyl chain not found"

    # Check for double bond between alpha and beta carbons
    bond = mol.GetBondBetweenAtoms(alpha_c.GetIdx(), beta_c.GetIdx())
    if bond is None or bond.GetBondType() != Chem.rdchem.BondType.DOUBLE:
        return False, "No double bond between positions 2 and 3 in acyl chain"

    return True, "2-enoyl-CoA identified with double bond at position 2 in acyl chain"