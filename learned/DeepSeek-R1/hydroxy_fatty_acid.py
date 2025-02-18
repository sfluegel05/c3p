"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: CHEBI hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmarts, GetAdjacencyMatrix

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is a fatty acid (carboxylic acid with aliphatic chain)
    carrying one or more hydroxy substituents on the aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find carboxylic acid groups (including deprotonated)
    carboxylic_acid = MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxylic_matches:
        return False, "No carboxylic acid group"

    # Check each carboxylic acid group's aliphatic chain
    for match in carboxylic_matches:
        cooh_carbon = match[0]
        visited = set()
        queue = [cooh_carbon]
        chain_carbons = set()

        # BFS to collect all aliphatic carbons connected to COOH carbon
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
            visited.add(current)
            atom = mol.GetAtomWithIdx(current)
            if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
                continue
            chain_carbons.add(current)
            # Traverse all non-aromatic carbon neighbors
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if neighbor.GetAtomicNum() != 6:
                    continue
                bond = mol.GetBondBetweenAtoms(current, n_idx)
                if bond.GetBondType() == Chem.BondType.AROMATIC:
                    continue
                if n_idx not in visited:
                    queue.append(n_idx)

        # Check chain length requirement (at least 4 carbons including COOH)
        if len(chain_carbons) < 4:
            continue

        # Check for hydroxyl groups on chain carbons (excluding COOH carbon)
        for carbon in chain_carbons:
            if carbon == cooh_carbon:
                continue
            atom = mol.GetAtomWithIdx(carbon)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    # Check if hydroxyl group is present
                    if neighbor.GetTotalNumHs() == 1 and \
                       mol.GetBondBetweenAtoms(carbon, neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                        return True, "Hydroxyl group on aliphatic chain of fatty acid"

    return False, "No hydroxyl group on aliphatic chain of fatty acid"