"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol is any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH,
    in which the carbon skeleton is composed of one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for alcohol group (-OH)
    OH_pattern = Chem.MolFromSmarts('[OX2H]')
    OH_matches = mol.GetSubstructMatches(OH_pattern)
    if not OH_matches:
        return False, "No hydroxyl group (-OH) found"

    # Find all carbon chains starting from an alcohol group
    for oh_match in OH_matches:
        oh_atom_idx = oh_match[0]
        oh_atom = mol.GetAtomWithIdx(oh_atom_idx)

        # Look for carbon atom attached to the hydroxyl oxygen
        neighbors = [n for n in oh_atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if not neighbors:
            continue  # No carbon attached to hydroxyl group
        start_atom = neighbors[0]

        # Traverse the longest carbon chain starting from the alcohol carbon
        path = Chem.GetLongestCarbonChain(mol, start_atom.GetIdx())
        if not path:
            continue  # Couldn't find a valid carbon chain

        # Check if the chain length corresponds to one or more isoprene units
        chain_length = len(path)
        if chain_length < 5 or chain_length % 5 != 0:
            continue  # Chain length does not correspond to complete isoprene units

        n_units = chain_length // 5  # Number of isoprene units

        # Check for methyl branches at correct positions
        is_prenol = True
        for i in range(n_units):
            unit_start_idx = path[i*5]
            unit_carbon = mol.GetAtomWithIdx(unit_start_idx)

            # The second carbon in the isoprene unit (position i*5 + 1)
            methyl_carbon_idx = path[i*5 + 1]
            methyl_carbon = mol.GetAtomWithIdx(methyl_carbon_idx)

            # Check for methyl branch attached to this carbon
            methyl_found = False
            for neighbor in methyl_carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in path:
                    # Found a carbon branch not in the main chain
                    if neighbor.GetDegree() == 1:
                        methyl_found = True
                        break
            if not methyl_found:
                is_prenol = False
                break  # Missing methyl group at expected position

        if is_prenol:
            return True, f"Molecule is a prenol with {n_units} isoprene unit(s)"

    return False, "Molecule does not conform to prenol structure"

# Helper function to get the longest carbon chain starting from a given atom
def GetLongestCarbonChain(mol, start_atom_idx):
    """
    Finds the longest carbon chain in the molecule starting from the given atom index.

    Args:
        mol: RDKit molecule object
        start_atom_idx: Index of the starting atom

    Returns:
        list: List of atom indices representing the longest carbon chain
    """
    def dfs(atom_idx, visited):
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        max_path = [atom_idx]
        for neighbor in atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx in visited:
                continue
            if neighbor.GetAtomicNum() != 6:
                continue  # Only consider carbon atoms
            path = dfs(nbr_idx, visited.copy())
            if len(path) > len(max_path):
                max_path = [atom_idx] + path
        return max_path

    return dfs(start_atom_idx, set())