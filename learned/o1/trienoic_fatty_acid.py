"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: trienoic fatty acid
Definition: Any polyunsaturated fatty acid that contains three double bonds.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a fatty acid with a terminal carboxylic acid group and
    exactly three carbon-carbon double bonds in its main unbranched aliphatic chain.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    carbox_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    if not carbox_matches:
        return False, "No carboxylic acid group found"

    # Identify terminal carboxylic acid group
    terminal_carboxyl_c = None
    for match in carbox_matches:
        c_idx = match[0]  # Carbonyl carbon atom index
        c_atom = mol.GetAtomWithIdx(c_idx)
        neighbors = c_atom.GetNeighbors()
        # Check if the carbon has only one neighbor besides the carboxyl oxygens
        non_oxygens = [n for n in neighbors if n.GetAtomicNum() != 8]
        if len(non_oxygens) == 1:
            terminal_carboxyl_c = c_idx
            break

    if terminal_carboxyl_c is None:
        return False, "No terminal carboxylic acid group found"

    # Use RDKit to find the longest carbon chain starting from the alpha carbon
    paths = Chem.rdmolops.GetShortestPathMatrix(mol)
    num_carbons = mol.GetNumAtoms(onlyExplicit=True)
    max_chain_length = 0
    main_chain = []

    # Identify all carbon chains starting from the alpha carbon
    alpha_carbon = None
    carboxyl_carbon = mol.GetAtomWithIdx(terminal_carboxyl_c)
    for neighbor in carboxyl_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != terminal_carboxyl_c:
            alpha_carbon = neighbor.GetIdx()
            break

    if alpha_carbon is None:
        return False, "No alpha carbon found adjacent to carboxylic group"

    # Use Depth-First Search to find the longest carbon chain
    def dfs(current_idx, visited, current_chain):
        nonlocal max_chain_length, main_chain
        visited.add(current_idx)
        current_chain.append(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        neighbors = [n.GetIdx() for n in current_atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() not in visited]
        is_terminal = len(neighbors) == 0
        if is_terminal:
            if len(current_chain) > max_chain_length:
                max_chain_length = len(current_chain)
                main_chain = current_chain.copy()
        for neighbor_idx in neighbors:
            # Allow oxygen substituents (e.g., hydroxy groups)
            neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
            dfs(neighbor_idx, visited.copy(), current_chain.copy())

    dfs(alpha_carbon, set(), [])

    if max_chain_length < 10:
        return False, f"Chain contains {max_chain_length + 1} carbon atoms, which is too short for a typical fatty acid"

    # Count the number of carbon-carbon double bonds in the main chain
    double_bonds = 0
    for i in range(len(main_chain) - 1):
        bond = mol.GetBondBetweenAtoms(main_chain[i], main_chain[i+1])
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bonds += 1

    if double_bonds != 3:
        return False, f"Contains {double_bonds} carbon-carbon double bonds in the main chain, needs exactly 3"

    # Check for branching (side chains attached to the main chain)
    for atom_idx in main_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        # Exclude main chain neighbors
        side_chains = [n for n in neighbors if n not in main_chain]
        if side_chains:
            return False, "Chain is branched, fatty acids are typically unbranched"

    # Verify the chain ends with the terminal carbon (no cycles)
    terminal_carbon = main_chain[-1]
    if mol.GetAtomWithIdx(terminal_carbon).IsInRing():
        return False, "Main chain forms a ring, fatty acids are acyclic"

    return True, "Molecule is a trienoic fatty acid with a terminal carboxylic acid group and exactly three C=C double bonds in the main unbranched aliphatic chain"