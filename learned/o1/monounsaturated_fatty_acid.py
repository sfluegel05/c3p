"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    A monounsaturated fatty acid has a carboxylic acid group and a long carbon chain
    with exactly one double or triple bond, and singly bonded carbon atoms in the rest of the chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (including deprotonated forms)
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O-,OH]")
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if not matches:
        return False, "No carboxylic acid group found"

    # Assume the first match is the carboxylic acid group
    carboxy_c_idx = matches[0][0]
    carboxy_o_idx = matches[0][1]

    # Find the carbon atom attached to the carboxylic acid carbon
    carboxy_c_atom = mol.GetAtomWithIdx(carboxy_c_idx)
    chain_c_atoms = set()

    # Use BFS to find the carbon chain connected to the carboxylic acid carbon
    visited = set()
    queue = [(carboxy_c_atom.GetIdx(), None)]
    while queue:
        current_idx, prev_idx = queue.pop(0)
        if current_idx in visited:
            continue
        visited.add(current_idx)
        atom = mol.GetAtomWithIdx(current_idx)
        if atom.GetAtomicNum() == 6:
            chain_c_atoms.add(current_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx != prev_idx and neighbor.GetAtomicNum() == 6:
                    queue.append((neighbor_idx, current_idx))
        else:
            continue

    # Exclude the carboxylic acid carbon
    chain_c_atoms.discard(carboxy_c_idx)
    if not chain_c_atoms:
        return False, "No carbon chain found attached to carboxylic acid"

    # Extract the subgraph of the chain
    chain_atoms = list(chain_c_atoms)
    chain_bonds = []
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx in chain_c_atoms and end_idx in chain_c_atoms:
            chain_bonds.append(bond)

    # Create a molecule of the chain
    chain_mol = Chem.PathToSubmol(mol, chain_bonds)

    # Count the number of double and triple bonds in the chain
    unsaturation_count = 0
    for bond in chain_mol.GetBonds():
        if bond.IsInRing():
            continue  # Exclude rings
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE or bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            unsaturation_count += 1

    # Include double bonds in rings (e.g., cyclopropene fatty acids)
    ring_bonds = chain_mol.GetRingInfo().BondRings()
    for ring in ring_bonds:
        for bond_idx in ring:
            bond = chain_mol.GetBondWithIdx(bond_idx)
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE or bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
                unsaturation_count += 1

    if unsaturation_count == 0:
        return False, "No double or triple bonds found in the chain"
    elif unsaturation_count > 1:
        return False, f"More than one double or triple bond found in the chain ({unsaturation_count} unsaturations)"

    # Check for additional functional groups
    allowed_elements = {1, 6, 8}  # H, C, O
    for atom in chain_mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_elements:
            return False, "Chain contains atoms other than carbon, hydrogen, and oxygen"
        if atom.GetAtomicNum() == 8:
            # Oxygen atoms should only be present in hydroxyl groups on the chain
            num_O_neighbors = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6)
            if num_O_neighbors > 1:
                return False, "Chain contains disallowed oxygen-containing functional groups"

    # Check chain length (optional threshold, e.g., chain length >= 4)
    if len(chain_c_atoms) < 4:
        return False, "Chain is too short to be a fatty acid"

    return True, "Molecule is a monounsaturated fatty acid"