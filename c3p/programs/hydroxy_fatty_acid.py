"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: hydroxy fatty acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is a fatty acid carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"

    # Get the carboxylic acid carbon atom indices
    carboxy_carbons = [match[0] for match in carboxy_matches]

    # Initialize variables to track the longest aliphatic chain
    max_chain_length = 0
    max_chain = []
    hydroxy_on_chain = False

    # Define patterns to exclude molecules with rings or non-aliphatic structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Molecule contains ring structures, not a typical fatty acid"

    # Atom indices of oxygen atoms in carboxylic acid groups
    carboxy_oxygen_indices = [match[1:] for match in carboxy_matches]
    carboxy_oxygen_indices = set([idx for sublist in carboxy_oxygen_indices for idx in sublist])

    # Iterate over each carboxylic acid group to find the longest aliphatic chain
    for carboxy_carbon_idx in carboxy_carbons:
        visited = set()
        stack = [(carboxy_carbon_idx, [carboxy_carbon_idx])]

        while stack:
            current_atom_idx, path = stack.pop()
            current_atom = mol.GetAtomWithIdx(current_atom_idx)

            if current_atom_idx in visited:
                continue
            visited.add(current_atom_idx)

            # Traverse to neighbor atoms
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in path:
                    continue  # Avoid cycles

                # Exclude carboxylic acid oxygens
                if neighbor_idx in carboxy_oxygen_indices:
                    continue

                bond = mol.GetBondBetweenAtoms(current_atom_idx, neighbor_idx)
                bond_type = bond.GetBondType()

                # Exclude non-single or double bonds (e.g., triple bonds)
                if bond_type not in [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE]:
                    continue

                neighbor_atomic_num = neighbor.GetAtomicNum()

                # Only consider carbon atoms in the chain
                if neighbor_atomic_num == 6:
                    new_path = path + [neighbor_idx]
                    stack.append((neighbor_idx, new_path))

                    # Check if this is the longest chain
                    if len(new_path) > max_chain_length:
                        max_chain_length = len(new_path)
                        max_chain = new_path

                # Check for hydroxy groups (-OH) attached to the chain carbons
                if neighbor_atomic_num == 8:
                    if neighbor.GetTotalDegree() == 1 and neighbor.GetTotalNumHs() >= 1:
                        bond_to_neighbor = mol.GetBondBetweenAtoms(current_atom_idx, neighbor_idx)
                        if bond_to_neighbor.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            hydroxy_on_chain = True

    if max_chain_length == 0:
        return False, "No aliphatic chain found connected to carboxylic acid group"

    if not hydroxy_on_chain:
        return False, "No hydroxy substituents found on the aliphatic chain"

    # Check for presence of other functional groups
    functional_groups = rdMolDescriptors.CalcNumHeteroatoms(mol)
    # Subtract the heteroatoms from carboxylic acid group and hydroxy substituents
    num_hydroxy_groups = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2H][CX4]')))
    num_carboxy_oxygens = len(carboxy_oxygen_indices)
    additional_heteroatoms = functional_groups - num_hydroxy_groups - num_carboxy_oxygens
    if additional_heteroatoms > 0:
        return False, "Molecule contains additional functional groups, not a typical fatty acid"

    # All checks passed
    return True, "Molecule is a hydroxy fatty acid"