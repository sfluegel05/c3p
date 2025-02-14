"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: unsaturated fatty acid
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid is any fatty acid containing at least one C=C or C#C bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid groups
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[OX1H1]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"

    # Identify all carbon-carbon double and triple bonds
    unsaturated_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() in (Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE):
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                unsaturated_bonds.append(bond)

    if not unsaturated_bonds:
        return False, "No carbon-carbon double or triple bonds in molecule"

    # Get indices of atoms involved in unsaturated bonds
    unsaturated_atom_indices = set()
    for bond in unsaturated_bonds:
        unsaturated_atom_indices.add(bond.GetBeginAtomIdx())
        unsaturated_atom_indices.add(bond.GetEndAtomIdx())

    # For each carboxyl carbon, check if there is a path to any unsaturated carbon atom
    for match in carboxy_matches:
        carboxyl_carbon_idx = match[0]
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

        # Use a BFS to find path to unsaturated carbon atoms
        visited = set()
        queue = [carboxyl_carbon_idx]

        while queue:
            current_idx = queue.pop(0)
            if current_idx in visited:
                continue
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            if current_idx in unsaturated_atom_indices and current_idx != carboxyl_carbon_idx:
                return True, "Molecule is an unsaturated fatty acid"

            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                neighbor_idx = neighbor.GetIdx()
                # Only consider carbon atoms
                if neighbor.GetAtomicNum() == 6 and neighbor_idx not in visited:
                    queue.append(neighbor_idx)

    return False, "Does not meet criteria for unsaturated fatty acid"