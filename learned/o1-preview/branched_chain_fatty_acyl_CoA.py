"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA
"""

from rdkit import Chem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA is a fatty acyl-CoA resulting from the condensation 
    of the thiol group of coenzyme A with the carboxy group of any branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    from rdkit import Chem
    from collections import deque

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find the thioester bond pattern (C(=O)-S-C)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Assume the first match is the thioester linkage
    carbonyl_c_idx = thioester_matches[0][0]
    sulfur_idx = thioester_matches[0][2]

    # Initialize traversal starting from the carbonyl carbon
    visited = set()
    acyl_chain_atoms = set()
    queue = deque()
    queue.append(carbonyl_c_idx)
    visited.add(sulfur_idx)  # Do not traverse into CoA moiety (sulfur side)

    while queue:
        atom_idx = queue.popleft()
        atom = mol.GetAtomWithIdx(atom_idx)
        acyl_chain_atoms.add(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited:
                bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
                # Avoid traversing into the CoA moiety via the sulfur atom
                if neighbor_idx != sulfur_idx:
                    visited.add(neighbor_idx)
                    queue.append(neighbor_idx)

    # Create a sub-molecule of the acyl chain
    acyl_chain_mol = Chem.PathToSubmol(mol, acyl_chain_atoms)

    # Check for rings in the acyl chain
    if acyl_chain_mol.GetRingInfo().NumRings() > 0:
        return False, "Acyl chain contains ring structures"

    # Check for heteroatoms other than oxygen (allowing for keto and hydroxyl groups)
    for atom in acyl_chain_mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 8]:  # Allow only H, C, O
            return False, f"Acyl chain contains heteroatom {atom.GetSymbol()}"

    # Check for branching in the acyl chain
    branching_points = 0
    for atom in acyl_chain_mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            neighbor_carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if len(neighbor_carbons) > 2:
                branching_points += 1

    if branching_points == 0:
        return False, "Acyl chain is not branched"
    else:
        return True, f"Acyl chain is branched with {branching_points} branching point(s)"