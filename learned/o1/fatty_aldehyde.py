"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: fatty aldehyde (CHEBI:35581)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde formally arising from reduction of the carboxylic acid group
    of its corresponding fatty acid, having a carbonyl group at one end of the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define the aldehyde group pattern (terminal aldehyde)
    # Aldehyde carbon (C=O) with single bond to hydrogen (implicit in SMILES) and single bond to carbon chain
    aldehyde_pattern = Chem.MolFromSmarts("[#6][CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

    if not aldehyde_matches:
        return False, "No aldehyde group found."

    # Check if aldehyde group is at the terminal position
    is_terminal_aldehyde = False
    for match in aldehyde_matches:
        aldehyde_c_idx = match[1]  # Index of aldehyde carbon
        aldehyde_c = mol.GetAtomWithIdx(aldehyde_c_idx)
        neighbors = aldehyde_c.GetNeighbors()

        # Aldehyde carbon should be connected to one carbon (chain) and one oxygen (double bond)
        if aldehyde_c.GetDegree() == 2:
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 6:
                    chain_atom_idx = neighbor.GetIdx()
                    is_terminal_aldehyde = True
                    break
        if is_terminal_aldehyde:
            break

    if not is_terminal_aldehyde:
        return False, "Aldehyde group is not at the terminal position."

    # Traverse the carbon chain starting from the carbon next to the aldehyde group
    visited = set()
    chain_length = 0

    def traverse_chain(atom_idx):
        nonlocal chain_length
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_length += 1
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx not in visited and neighbor.GetAtomicNum() in [6, 1]:
                    traverse_chain(n_idx)

    traverse_chain(chain_atom_idx)

    # A minimal fatty aldehyde should have at least 4 carbons in the chain
    if chain_length < 4:
        return False, f"Carbon chain is too short ({chain_length} carbons)."

    # Check for the presence of disallowed elements (only common bioelements allowed)
    allowed_atomic_nums = {1, 6, 7, 8, 15, 16, 17}  # H, C, N, O, P, S, Cl
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element with atomic number {atom.GetAtomicNum()}."

    return True, "Molecule is a fatty aldehyde with a terminal aldehyde group and appropriate carbon chain."