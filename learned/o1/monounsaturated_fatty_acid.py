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
    Hydroxyl groups may be present on the chain.

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

    # Exclude carboxyl oxygen atoms
    carboxy_o_indices = matches[0][1:]

    # Use BFS to find the carbon chain connected to the carboxylic acid carbon
    chain_atom_indices = set()
    visited = set()
    queue = [carboxy_c_idx]
    while queue:
        current_idx = queue.pop(0)
        if current_idx in visited:
            continue
        visited.add(current_idx)
        atom = mol.GetAtomWithIdx(current_idx)

        # Include carbon and oxygen atoms (allowing for hydroxyl groups)
        if atom.GetAtomicNum() in (6, 8):
            chain_atom_indices.add(current_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Exclude carboxyl oxygen atoms
                if neighbor_idx not in visited and neighbor_idx not in carboxy_o_indices:
                    queue.append(neighbor_idx)
        else:
            # Stop extending the chain at atoms other than carbon or oxygen
            continue

    # Remove the carboxylic acid carbon and oxygens from the chain atoms
    chain_atom_indices.discard(carboxy_c_idx)
    for idx in carboxy_o_indices:
        chain_atom_indices.discard(idx)

    if not chain_atom_indices:
        return False, "No carbon chain found attached to carboxylic acid"

    # Collect bonds in the chain
    chain_bond_indices = []
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx in chain_atom_indices and end_idx in chain_atom_indices:
            chain_bond_indices.append(bond.GetIdx())

    # Create a molecule of the chain
    chain_mol = Chem.PathToSubmol(mol, chain_bond_indices)

    # Count the number of double and triple bonds in the chain
    unsaturation_count = 0
    for bond in chain_mol.GetBonds():
        bond_type = bond.GetBondType()
        if bond_type == Chem.rdchem.BondType.DOUBLE or bond_type == Chem.rdchem.BondType.TRIPLE:
            unsaturation_count += 1

    if unsaturation_count == 0:
        return False, "No double or triple bonds found in the chain"
    elif unsaturation_count > 1:
        return False, f"More than one double or triple bond found in the chain ({unsaturation_count} unsaturations)"

    # Check that the rest of the chain contains singly bonded carbons (allowing for hydroxyl groups)
    for atom in chain_mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:  # Carbon atom
            # Exclude the atom if it is part of the unsaturation
            num_double_bonds = sum(1 for bond in atom.GetBonds() if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE])
            if num_double_bonds > 1:
                return False, "Carbon atom in chain has more than one unsaturation"
        elif atomic_num == 8:  # Oxygen atom
            # Allow oxygen atoms in hydroxyl groups (attached to carbon)
            num_oxygen_bonds = len(atom.GetBonds())
            if num_oxygen_bonds > 1:
                return False, "Oxygen atom in chain with invalid bonding"
        else:
            return False, "Chain contains atoms other than carbon and oxygen"

    # Check chain length (e.g., chain length >= 4 carbons)
    num_carbons_in_chain = sum(1 for idx in chain_atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if num_carbons_in_chain < 4:
        return False, "Chain is too short to be a fatty acid"

    return True, "Molecule is a monounsaturated fatty acid"