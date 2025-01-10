"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for phosphate ester group attached via oxygen
    phosphate_pattern = Chem.MolFromSmarts("[O]-P(=O)([O])[O]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)

    if not phosphate_matches:
        return False, "No phosphate ester group found"

    # For each phosphate ester group, check if it is attached to a monosaccharide unit
    for match in phosphate_matches:
        # match[0] is the index of the oxygen atom connected to P
        oxygen_idx = match[0]
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)

        # Get the atom connected to the oxygen (excluding phosphorus)
        neighbors = [atom for atom in oxygen_atom.GetNeighbors() if atom.GetAtomicNum() != 15]
        if not neighbors:
            continue
        carbon_atom = neighbors[0]
        if carbon_atom.GetAtomicNum() != 6:
            continue  # Not connected to carbon

        # Now, starting from the carbon atom, collect the connected atoms up to a depth limit
        visited = set()
        atoms_to_visit = [carbon_atom.GetIdx()]
        subgraph_atoms = set()

        while atoms_to_visit:
            current_idx = atoms_to_visit.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)
            subgraph_atoms.add(current_idx)

            atomic_num = current_atom.GetAtomicNum()
            if atomic_num not in [6, 8]:  # Only consider C and O atoms
                continue  # Stop traversal at heteroatoms

            # Add neighbor atoms to visit
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                neighbor_atomic_num = neighbor.GetAtomicNum()
                if neighbor_idx not in visited and neighbor_atomic_num in [6, 8]:
                    atoms_to_visit.append(neighbor_idx)

        # Analyze the subgraph to determine if it's a monosaccharide
        num_carbons = 0
        num_oxygens = 0
        for idx in subgraph_atoms:
            atom = mol.GetAtomWithIdx(idx)
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 6:
                num_carbons += 1
            elif atomic_num == 8:
                num_oxygens +=1

        # Monosaccharides have 3 to 7 carbons and multiple hydroxyl groups
        if 3 <= num_carbons <=7:
            if num_oxygens >= num_carbons - 1:
                return True, "Contains monosaccharide unit with phosphate ester group attached to hydroxyl group"

    return False, "Phosphate ester group not attached to monosaccharide unit"

__metadata__ = {
    'chemical_class': {
        'name': 'phospho sugar',
        'definition': 'Any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.'
    }
}