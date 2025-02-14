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
    A monounsaturated fatty acid has a carboxylic acid group and a long unbranched carbon chain
    with exactly one double or triple bond (excluding the carboxylic acid group), and singly bonded carbon atoms in the rest of the chain.
    Hydroxyl groups may be present on the chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts("C(=O)O")
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if not matches:
        return False, "No carboxylic acid group found"

    if len(matches) > 1:
        return False, "More than one carboxylic acid group found"

    # Get the carboxylic acid carbon atom index
    carboxy_c_idx = matches[0][0]

    # Recursively find the longest unbranched carbon chain starting from carboxylic carbon
    def find_longest_chain(current_idx, visited):
        atom = mol.GetAtomWithIdx(current_idx)
        if atom.GetAtomicNum() != 6:
            return [], []

        visited.add(current_idx)
        max_chain = [current_idx]
        max_bonds = []
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(current_idx, neighbor_idx)
            if neighbor_idx not in visited:
                neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
                if neighbor_atom.GetAtomicNum() == 6 or (neighbor_atom.GetAtomicNum() == 8 and neighbor_atom.GetDegree() == 1):
                    # Exclude branching in carbon atoms (degree > 2)
                    if neighbor_atom.GetAtomicNum() == 6 and neighbor_atom.GetDegree() > 2:
                        continue
                    chain, bonds = find_longest_chain(neighbor_idx, visited)
                    chain = [current_idx] + chain
                    bonds = [bond.GetIdx()] + bonds
                    if len(chain) > len(max_chain):
                        max_chain = chain
                        max_bonds = bonds
        visited.remove(current_idx)
        return max_chain, max_bonds

    longest_chain, chain_bonds = find_longest_chain(carboxy_c_idx, set())
    if len(longest_chain) <= 1:
        return False, "No carbon chain found attached to carboxylic acid"

    chain_mol = Chem.PathToSubmol(mol, chain_bonds)

    # Check for rings in the chain
    sssr = Chem.GetSSSR(chain_mol)
    if sssr > 0:
        return False, "Chain contains rings"

    # Count the number of double and triple bonds in the chain, excluding the carboxylic acid
    unsaturation_count = 0
    for bond in chain_mol.GetBonds():
        bond_type = bond.GetBondType()
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if carboxy_c_idx in [begin_idx, end_idx]:
            continue  # Skip bonds involving carboxylic acid carbon
        if bond_type in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            unsaturation_count += 1

    if unsaturation_count == 0:
        return False, "No double or triple bonds found in the chain"
    elif unsaturation_count > 1:
        return False, f"More than one double or triple bond found in the chain ({unsaturation_count} unsaturations)"

    # Check that the rest of the chain contains singly bonded carbons (allowing for hydroxyl groups)
    for atom in chain_mol.GetAtoms():
        idx = atom.GetIdx()
        atomic_num = atom.GetAtomicNum()
        if idx == carboxy_c_idx:
            continue  # Skip carboxylic acid carbon
        if atomic_num == 6:  # Carbon atom
            num_multiple_bonds = sum(1 for bond in atom.GetBonds()
                                     if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]
                                     and bond.GetOtherAtom(atom).GetIdx() in longest_chain)
            if num_multiple_bonds > 1:
                return False, "Carbon atom in chain has more than one unsaturation"
            if atom.GetDegree() > 2:
                return False, "Chain is branched at carbon atom"
        elif atomic_num == 8:  # Oxygen atom
            # Allow hydroxyl groups (oxygen with degree 1)
            if atom.GetDegree() != 1:
                return False, "Oxygen atom in chain with invalid bonding"
        else:
            return False, "Chain contains atoms other than carbon and oxygen"

    # Optionally, adjust chain length criteria
    num_carbons_in_chain = sum(1 for idx in longest_chain if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if num_carbons_in_chain < 2:
        return False, "Chain is too short to be a fatty acid"

    return True, "Molecule is a monounsaturated fatty acid"