"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: ultra-long-chain fatty acid
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    An ultra-long-chain fatty acid is defined as any very long-chain fatty acid which has a chain length greater than C27.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for exactly one carboxylic acid functional group [CX3](=O)[OX1H1]
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX1H1]")
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(matches) != 1:
        return False, f"Expected one carboxylic acid group, found {len(matches)}"

    # Get the carboxyl carbon atom index
    carboxyl_carbon_idx = matches[0][0]

    # Identify all terminal carbon atoms (degree 1 carbons)
    terminal_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and len(atom.GetNeighbors()) == 1:
            terminal_carbons.append(atom.GetIdx())

    # Check if there are terminal carbons
    if not terminal_carbons:
        return False, "No terminal carbon atoms found"

    # Function to check if a path consists only of carbons
    def path_is_linear_carbon_chain(path):
        bond_types = [mol.GetBondBetweenAtoms(path[i], path[i+1]).GetBondType() for i in range(len(path)-1)]
        atoms = [mol.GetAtomWithIdx(idx) for idx in path]
        # Check if all atoms are carbons and all bonds are single
        is_linear = all(atom.GetAtomicNum() == 6 for atom in atoms) and all(bond == Chem.rdchem.BondType.SINGLE for bond in bond_types)
        return is_linear

    # Find the longest linear carbon chain from carboxyl carbon to terminal carbons
    from rdkit.Chem import rdmolops

    max_chain_length = 0
    for term_idx in terminal_carbons:
        # Get all paths between carboxyl carbon and terminal carbon
        paths = rdmolops.GetAllSimplePaths(mol, carboxyl_carbon_idx, term_idx)
        for path in paths:
            if path_is_linear_carbon_chain(path):
                chain_length = len(path)
                if chain_length > max_chain_length:
                    max_chain_length = chain_length

    if max_chain_length == 0:
        return False, "No linear carbon chain found from carboxyl carbon to terminal carbon"

    # Subtract 1 to get the number of bonds (chain length in carbons)
    chain_length = max_chain_length

    if chain_length > 27:
        return True, f"Longest linear carbon chain is {chain_length} carbons"
    else:
        return False, f"Longest linear carbon chain is {chain_length} carbons, which is not greater than 27"