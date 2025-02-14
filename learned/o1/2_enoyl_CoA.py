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

    # Define a general SMARTS pattern for Coenzyme A (CoA) moiety
    # CoA contains a phosphoadenosine diphosphate (ADP) and pantetheine moiety
    coa_smarts = '[#8]-P(=O)([O-])-[O]-[#6]-[#5]1[#6](=[#7])[#7]=[#6]([#7]=[#6]1)-[#8]-[#5]2[#6](=[#7])[#7]=[#6]([#7]=[#6]2)-[#7]'
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A (CoA) moiety not found"

    # Define SMARTS pattern for thioester linkage: C(=O)S
    thioester_smarts = 'C(=O)S'
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Assume only one thioester linkage (acyl chain attached to CoA)
    thioester_match = thioester_matches[0]
    carbonyl_c_idx = thioester_match[0]  # Carbonyl carbon index
    sulfur_idx = thioester_match[2]      # Sulfur atom index

    # Traverse the acyl chain starting from the carbonyl carbon
    # Collect acyl chain atoms until reaching the sulfur atom
    acyl_chain_atoms = []
    visited = set()
    stack = [carbonyl_c_idx]
    while stack:
        atom_idx = stack.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            continue  # Only consider carbon atoms in the acyl chain
        acyl_chain_atoms.append(atom_idx)
        # Add neighboring carbons to stack
        for neighbor in atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 6 and nbr_idx not in visited:
                stack.append(nbr_idx)

    # Sort acyl chain atoms based on their positions in the molecule
    acyl_chain_atoms.sort()

    # Identify the positions of the acyl chain atoms relative to the carbonyl carbon
    atom_positions = {carbonyl_c_idx: 1}
    for idx in acyl_chain_atoms:
        if idx == carbonyl_c_idx:
            continue
        try:
            path_length = Chem.GetShortestPath(mol, carbonyl_c_idx, idx)
            atom_positions[idx] = len(path_length)
        except:
            continue

    # Check for double bond between positions 2 and 3
    double_bond_found = False
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        bond_order = bond.GetBondTypeAsDouble()
        if bond_order == 2.0:
            pos1 = atom_positions.get(begin_idx, None)
            pos2 = atom_positions.get(end_idx, None)
            if pos1 in [2,3] and pos2 in [2,3] and pos1 != pos2:
                double_bond_found = True
                break

    if not double_bond_found:
        return False, "No double bond between positions 2 and 3 in the acyl chain"

    return True, "Molecule is a 2-enoyl-CoA with a double bond between positions 2 and 3"