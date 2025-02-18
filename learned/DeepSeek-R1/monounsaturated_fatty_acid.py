"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid.
    A monounsaturated fatty acid has a carboxylic acid group and exactly one double or triple bond in the carbon chain.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for carboxylic acid group
    carboxyl = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl):
        return False, "No carboxylic acid group"

    # Get the main chain (longest carbon chain starting from the carboxylic acid)
    matches = mol.GetSubstructMatches(carboxyl)
    if not matches:
        return False, "No carboxylic acid"
    # Assuming the first match is the COOH group; get the adjacent carbon
    cooh_carbon = matches[0][0]
    chain = []
    visited = set()
    stack = [(cooh_carbon, 0)]  # (atom index, chain length)
    max_length = 0
    while stack:
        atom_idx, length = stack.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:  # Focus on carbon chain
            continue
        chain.append(atom_idx)
        neighbors = atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                stack.append((neighbor.GetIdx(), length + 1))
        if length > max_length:
            max_length = length

    # Count double/triple bonds in the chain
    unsat_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in chain and bond.GetEndAtomIdx() in chain:
            if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
                unsat_bonds += 1

    if unsat_bonds != 1:
        return False, f"Found {unsat_bonds} unsaturated bonds, need exactly 1"

    # Check chain length (at least 4 carbons excluding COOH)
    if max_length < 4:
        return False, "Chain too short"

    # Check all other bonds in the chain are single
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in chain and bond.GetEndAtomIdx() in chain:
            if bond.GetBondType() not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
                return False, "Invalid bond type in chain"
            # Ensure only one unsaturation
            if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
                if unsat_bonds != 1:
                    return False, "Multiple unsaturated bonds"

    return True, "Monounsaturated fatty acid with one unsaturated bond and carboxylic acid group"