"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is a fatty acid which has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid group (C(=O)O[H] or deprotonated form)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H0-1]')
    if carboxylic_acid_pattern is None:
        return False, "Invalid SMARTS pattern for carboxylic acid"
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Exclude molecules with ester or amide functional groups connected to the chain
    ester_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2][!H]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if ester_matches:
        return False, "Ester functional group(s) found"

    amide_pattern = Chem.MolFromSmarts('[CX3](=O)[NX3][!$([NX3][H2])]')
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if amide_matches:
        return False, "Amide functional group(s) found"

    # Identify the aliphatic chain attached to the carboxyl group
    # We will find the longest path starting from the carbonyl carbon
    fatty_acid = False
    max_chain_length = 0
    reason = ""

    for match in carboxylic_acid_matches:
        carbonyl_c_idx = match[0]
        carbonyl_c_atom = mol.GetAtomWithIdx(carbonyl_c_idx)

        # Perform a search for the longest aliphatic chain starting from the carbonyl carbon
        visited = set()
        stack = [(carbonyl_c_idx, 0)]  # (atom_idx, chain_length)

        while stack:
            current_idx, chain_length = stack.pop()
            visited.add(current_idx)
            current_atom = mol.GetAtomWithIdx(current_idx)

            # Exclude non-carbon atoms
            if current_atom.GetAtomicNum() != 6:
                continue

            # Update maximum chain length
            if chain_length > max_chain_length:
                max_chain_length = chain_length

            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    bond = mol.GetBondBetweenAtoms(current_idx, neighbor_idx)
                    # Exclude cycles
                    if bond.IsInRing():
                        continue
                    # Exclude backward to carboxyl oxygens
                    if neighbor.GetAtomicNum() == 8:
                        continue
                    # Include only single, double, or triple bonds (to account for unsaturation)
                    if bond.GetBondType() in (rdchem.BondType.SINGLE, rdchem.BondType.DOUBLE, rdchem.BondType.TRIPLE):
                        stack.append((neighbor_idx, chain_length + 1))

    if max_chain_length > 22:
        fatty_acid = True
        reason = f"Longest carbon chain length is {max_chain_length}, which is greater than 22"
    else:
        reason = f"Longest carbon chain length is {max_chain_length}, which is not greater than 22"

    return fatty_acid, reason