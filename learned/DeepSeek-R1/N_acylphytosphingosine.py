"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: CHEBI:78231 N-acylphytosphingosine
"""
from rdkit import Chem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide with a phytosphingosine backbone (at least three hydroxyl groups)
    and a fatty acyl group attached via an amide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find amide groups [N connected to carbonyl]
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"

    for amide_match in amide_matches:
        nitrogen_idx = amide_match[0]
        carbonyl_idx = amide_match[1]
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        carbonyl = mol.GetAtomWithIdx(carbonyl_idx)

        # Check sphingoid base for at least three hydroxyl groups
        visited = set([carbonyl_idx])
        stack = [nitrogen]
        hydroxyl_count = 0

        while stack:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())

            # Check for hydroxyl groups attached to this atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    hydroxyl_count += 1

            # Add neighbors except those leading to the acyl group
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited:
                    stack.append(neighbor)

        if hydroxyl_count < 3:
            continue  # Not enough hydroxyls in sphingoid base

        # Check acyl chain length (R in R-C(=O)-N...)
        acyl_chain_carbons = 0
        visited_acyl = set()
        # Get R groups attached to carbonyl (excluding nitrogen)
        r_groups = []
        for neighbor in carbonyl.GetNeighbors():
            if neighbor.GetIdx() != nitrogen_idx:
                r_groups.append(neighbor)

        # Traverse each R group to count carbons
        for r in r_groups:
            stack_acyl = [r]
            while stack_acyl:
                atom = stack_acyl.pop()
                if atom.GetIdx() in visited_acyl:
                    continue
                visited_acyl.add(atom.GetIdx())
                if atom.GetAtomicNum() == 6:
                    acyl_chain_carbons += 1
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited_acyl:
                        stack_acyl.append(neighbor)

        if acyl_chain_carbons >= 4:  # Minimum of 4 carbons for fatty acyl
            return True, f"Phytosphingosine backbone with {hydroxyl_count} hydroxyls and {acyl_chain_carbons}-carbon acyl chain"

    return False, "Does not meet N-acylphytosphingosine criteria"