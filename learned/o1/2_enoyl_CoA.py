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

    # Define SMARTS patterns for key features
    # Thioester linkage: C(=O)-S-C
    thioester_smarts = '[CX3](=O)[SX2][#6]'
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)

    # Adenine ring in CoA
    adenine_smarts = 'n1c([nH])cnc1N'
    adenine_pattern = Chem.MolFromSmarts(adenine_smarts)

    # Ribose sugar connected to adenine
    ribose_smarts = 'OC[C@H]1O[C@H](CO)[C@@H](O)[C@H]1O'
    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)

    # Diphosphate linkage
    diphosphate_smarts = 'OP(=O)(O)OP(=O)(O)O'
    diphosphate_pattern = Chem.MolFromSmarts(diphosphate_smarts)

    # Pantetheine moiety (portion of CoA connected to thioester)
    pantetheine_smarts = '[NX3][CX3](=O)[CX4][CX3](=O)[NX3][CX4][CX4][SX2]'
    pantetheine_pattern = Chem.MolFromSmarts(pantetheine_smarts)

    # Check for thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Assume the first thioester linkage is the one we're interested in
    thioester_match = thioester_matches[0]
    carbonyl_c_idx = thioester_match[0]  # Carbonyl carbon index
    sulfur_idx = thioester_match[2]      # Sulfur atom index

    # Check for CoA substructures connected to sulfur
    sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
    coa_found = False

    # Perform recursive search from sulfur atom to find CoA moiety
    atoms_to_visit = [sulfur_atom]
    visited_atoms = set()
    while atoms_to_visit:
        current_atom = atoms_to_visit.pop()
        if current_atom.GetIdx() in visited_atoms:
            continue
        visited_atoms.add(current_atom.GetIdx())

        # Check for adenine ring
        if current_atom.HasSubstructMatch(adenine_pattern):
            coa_found = True
            break

        # Add neighbors to visit
        for neighbor in current_atom.GetNeighbors():
            atoms_to_visit.append(neighbor)

    if not coa_found:
        return False, "Coenzyme A (CoA) moiety not found"

    # Identify the acyl chain connected to the carbonyl carbon
    # We will traverse away from the carbonyl carbon, avoiding the sulfur atom
    acyl_chain = []
    visited = set()
    stack = [carbonyl_c_idx]
    position_map = {carbonyl_c_idx: 1}  # Map atom idx to position in acyl chain

    while stack:
        atom_idx = stack.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        acyl_chain.append(atom_idx)
        current_position = position_map[atom_idx]
        for neighbor in atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx == sulfur_idx:
                continue  # Skip the sulfur atom leading to CoA
            if nbr_idx in visited:
                continue
            # Ensure we are moving along the acyl chain (carbons and hydrogens)
            nbr_atom = mol.GetAtomWithIdx(nbr_idx)
            if nbr_atom.GetAtomicNum() in [1]:  # Skip hydrogens
                continue
            # Assign position number
            position_map[nbr_idx] = current_position + 1
            stack.append(nbr_idx)

    # Check if acyl chain is long enough
    if len(acyl_chain) < 3:
        return False, "Acyl chain is too short"

    # Check for double bond between positions 2 and 3
    double_bond_found = False
    for bond in mol.GetBonds():
        bond_order = bond.GetBondTypeAsDouble()
        if bond_order == 2.0:
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            pos1 = position_map.get(begin_idx)
            pos2 = position_map.get(end_idx)
            if pos1 is None or pos2 is None:
                continue
            if {pos1, pos2} == {2, 3}:
                double_bond_found = True
                break

    if not double_bond_found:
        return False, "No double bond between positions 2 and 3 in the acyl chain"

    return True, "Molecule is a 2-enoyl-CoA with a double bond between positions 2 and 3"