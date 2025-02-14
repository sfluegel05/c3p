"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is an unsaturated fatty acyl-CoA with a double bond between positions 2 and 3
    of the acyl chain attached via a thioester linkage to Coenzyme A.

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

    # Define CoA molecule from its SMILES
    smiles_coa = 'CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C2N=CN=C3N)O)O)OC(=O)CCNC(=O)CCSC'
    coa_mol = Chem.MolFromSmiles(smiles_coa)
    if coa_mol is None:
        return False, "Failed to create CoA molecule"

    # Check for CoA substructure
    if not mol.HasSubstructMatch(coa_mol):
        return False, "Coenzyme A (CoA) moiety not found"

    # Define thioester linkage pattern: C(=O)S
    thioester_smarts = 'C(=O)S'
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Assume the first thioester linkage is the one we're interested in
    thioester_match = thioester_matches[0]
    carbonyl_c_idx = thioester_match[0]  # Carbonyl carbon index
    sulfur_idx = thioester_match[2]      # Sulfur atom index

    # Find the acyl chain connected to the carbonyl carbon
    # Traverse away from the carbonyl carbon, avoiding the sulfur atom
    acyl_chain = []
    visited = set()
    stack = [carbonyl_c_idx]
    while stack:
        atom_idx = stack.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        acyl_chain.append(atom_idx)
        for neighbor in atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx == sulfur_idx:
                continue  # Skip the sulfur atom leading to CoA
            if nbr_idx not in visited:
                stack.append(nbr_idx)

    # Remove the carbonyl carbon (position 1)
    acyl_chain.remove(carbonyl_c_idx)

    # Check that acyl chain has at least 2 carbons
    if len(acyl_chain) < 2:
        return False, "Acyl chain is too short"

    # Build a map of positions in the acyl chain
    atom_positions = {carbonyl_c_idx: 1}
    current_positions = [carbonyl_c_idx]
    position = 1
    while current_positions:
        next_positions = []
        position += 1
        for idx in current_positions:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if nbr_idx in atom_positions:
                    continue
                if nbr_idx == sulfur_idx:
                    continue  # Skip sulfur atom
                if nbr_idx in acyl_chain:
                    atom_positions[nbr_idx] = position
                    next_positions.append(nbr_idx)
        current_positions = next_positions

    # Check for double bond between positions 2 and 3
    double_bond_found = False
    for bond in mol.GetBonds():
        bond_order = bond.GetBondTypeAsDouble()
        if bond_order == 2.0:
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            pos1 = atom_positions.get(begin_idx)
            pos2 = atom_positions.get(end_idx)
            if {pos1, pos2} == {2, 3}:
                double_bond_found = True
                break

    if not double_bond_found:
        return False, "No double bond between positions 2 and 3 in the acyl chain"

    return True, "Molecule is a 2-enoyl-CoA with a double bond between positions 2 and 3"