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
    thioester_pattern = Chem.MolFromSmarts("C(=O)[SX1]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Identify adenine moiety (part of CoA)
    adenine_smarts = "n1cnc2c1ncnc2N"  # SMARTS pattern for adenine
    adenine_mol = Chem.MolFromSmarts(adenine_smarts)
    adenine_matches = mol.GetSubstructMatches(adenine_mol)
    if not adenine_matches:
        return False, "Adenine moiety not found (CoA not detected)"

    # Get indices of adenine atoms
    adenine_atom_indices = set()
    for match in adenine_matches:
        adenine_atom_indices.update(match)

    # For each thioester linkage, check if sulfur is connected to adenine
    for thioester_match in thioester_matches:
        carbonyl_c_idx = thioester_match[0]
        sulfur_idx = thioester_match[2]
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)

        # Check if sulfur is connected to adenine moiety
        connected_to_adenine = False
        for adenine_atom_idx in adenine_atom_indices:
            path = Chem.rdmolops.GetShortestPath(mol, sulfur_idx, adenine_atom_idx)
            if path:
                connected_to_adenine = True
                break
        if not connected_to_adenine:
            continue  # Try next thioester linkage

        # Identify alpha carbon (attached to carbonyl carbon but not sulfur)
        carbonyl_c_atom = mol.GetAtomWithIdx(carbonyl_c_idx)
        alpha_c = None
        for neighbor in carbonyl_c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != sulfur_idx:
                alpha_c = neighbor
                break
        if alpha_c is None:
            continue  # No alpha carbon found in this thioester linkage

        # Identify beta carbon (double-bonded to alpha carbon)
        beta_c = None
        for neighbor in alpha_c.GetNeighbors():
            if neighbor.GetIdx() == carbonyl_c_idx:
                continue
            if neighbor.GetAtomicNum() == 6:
                bond = mol.GetBondBetweenAtoms(alpha_c.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    beta_c = neighbor
                    break
        if beta_c is None:
            continue  # No beta carbon double-bonded to alpha carbon

        # Ensure double bond is not part of a ring
        if alpha_c.IsInRing() or beta_c.IsInRing():
            continue  # Double bond is part of a ring

        # Check that the acyl chain continues beyond beta carbon
        chain_continues = False
        for neighbor in beta_c.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != alpha_c.GetIdx():
                chain_continues = True
                break
        if not chain_continues:
            continue  # Acyl chain does not continue beyond beta carbon

        # If all checks pass, return True
        return True, "2-enoyl-CoA identified with double bond between positions 2 and 3 in acyl chain"

    return False, "Did not match criteria for 2-enoyl-CoA"