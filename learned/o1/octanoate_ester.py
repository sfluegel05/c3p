"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any ester where the acyl component is octanoic acid (caprylic acid).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester functional group pattern
    ester_pattern = Chem.MolFromSmarts("[$([CX3](=O)[OX2H0])]")  # Ester group pattern

    # Find all ester groups in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Iterate over each ester group
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Index of the carbonyl carbon
        ester_o_idx = match[2]     # Index of the ester oxygen

        # Initialize variables for traversal
        visited = set()
        stack = [(carbonyl_c_idx, None)]  # Stack of (current atom idx, previous atom idx)
        acyl_chain_length = 0
        branching = False

        # Traverse the acyl chain starting from the carbonyl carbon
        while stack:
            current_atom_idx, previous_atom_idx = stack.pop()
            if current_atom_idx in visited:
                continue
            visited.add(current_atom_idx)
            current_atom = mol.GetAtomWithIdx(current_atom_idx)

            # Count only carbon atoms (excluding the carbonyl carbon)
            if current_atom.GetAtomicNum() == 6 and current_atom_idx != carbonyl_c_idx:
                acyl_chain_length += 1

            # Get neighbors excluding the previous atom to avoid backtracking
            neighbors = [nbr for nbr in current_atom.GetNeighbors() if nbr.GetIdx() != previous_atom_idx]

            # Check for branching (should only have two neighbors in acyl chain)
            if len(neighbors) > 1 and current_atom_idx != carbonyl_c_idx:
                branching = True
                break

            for neighbor in neighbors:
                nbr_idx = neighbor.GetIdx()
                nbr_atom = mol.GetAtomWithIdx(nbr_idx)
                nbr_atom_num = nbr_atom.GetAtomicNum()

                # Exclude ester oxygen and any heteroatoms (to ensure unbranched carbon chain)
                if nbr_idx == ester_o_idx:
                    continue  # Skip ester oxygen
                elif nbr_atom_num == 6:
                    stack.append((nbr_idx, current_atom_idx))
                else:
                    branching = True
                    break

            if branching:
                break

        # Check if acyl chain meets the criteria for octanoic acid
        if branching:
            continue  # This ester has a branched acyl chain; check next ester
        if acyl_chain_length != 7:
            continue  # Acyl chain does not have 7 carbons (excluding carbonyl carbon)
        else:
            # Check for saturation (no double bonds in acyl chain)
            acyl_chain_atoms = list(visited)
            unsaturation = False
            for atom_idx in acyl_chain_atoms:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6:
                    for bond in atom.GetBonds():
                        if bond.GetBondType() == rdchem.BondType.DOUBLE and bond.GetOtherAtomIdx(atom_idx) in acyl_chain_atoms:
                            unsaturation = True
                            break
                if unsaturation:
                    break
            if unsaturation:
                continue  # Acyl chain is unsaturated; check next ester

            # Found an octanoate ester
            return True, "Molecule contains an octanoate ester group"

    # No octanoate ester groups found
    return False, "No octanoate ester groups found"