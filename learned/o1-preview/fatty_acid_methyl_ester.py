"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: CHEBI:XXXX fatty acid methyl ester
"""
from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is the carboxylic ester obtained by the formal condensation
    of a fatty acid with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the methyl ester SMARTS pattern: [C](=O)O[C]
    methyl_ester_smarts = "[CX3](=O)OC"
    methyl_ester_pattern = Chem.MolFromSmarts(methyl_ester_smarts)

    # Find all methyl ester groups in the molecule
    ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    if not ester_matches:
        return False, "No methyl ester group found"

    # For each methyl ester group found
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon
        ester_o_idx = match[1]     # Ester oxygen
        methyl_c_idx = match[2]    # Methyl carbon

        # Get the carbonyl carbon atom
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)

        # Get neighbors of the carbonyl carbon excluding the ester oxygen
        chain_atoms = []
        visited = set([carbonyl_c_idx, ester_o_idx, methyl_c_idx])
        stack = []

        for neighbor in carbonyl_c.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited:
                stack.append(neighbor)

        # Traverse the carbon chain attached to the carbonyl carbon
        while stack:
            atom = stack.pop()
            atom_idx = atom.GetIdx()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atomic_num = atom.GetAtomicNum()

            if atomic_num == 6:  # Carbon atom
                chain_atoms.append(atom)
                # Add neighboring atoms to the stack
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited:
                        stack.append(neighbor)
            elif atomic_num in [1, 8]:  # Hydrogen or Oxygen (allow hydroxyl groups)
                continue
            else:
                # Contains other heteroatoms (e.g., N, S), disqualify
                return False, "Chain contains heteroatoms, not a fatty acid chain"

        # Check if the carbon chain length is at least 4 (typical for fatty acids)
        chain_length = len(chain_atoms)
        if chain_length >= 4:
            return True, "Molecule is a fatty acid methyl ester"

    # If no suitable chain is found
    return False, "Methyl ester group not attached to a fatty acid chain of sufficient length"