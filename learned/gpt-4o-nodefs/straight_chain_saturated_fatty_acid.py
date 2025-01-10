"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    A straight-chain saturated fatty acid is characterized by a long continuous chain of carbons
    with a terminal carboxylic acid group, and no branching or cyclic structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the carboxylic acid group pattern at one end
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylic acid group found"

    # Ensure no double or triple bonds in the main carbon chain
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False, "Carbon chain contains double or triple bonds"

    # Track visited atoms and ensure a continuous linear path from carboxylic group
    stack = []
    visited = set()

    # Get carboxylic group carbon
    carboxyl_match = mol.GetSubstructMatch(carboxylate_pattern)
    if len(carboxyl_match) < 1:
        return False, "No carboxylic group found"
    start_atom = mol.GetAtomWithIdx(carboxyl_match[0])
    stack.append((start_atom, None))  # (current atom, previous atom)

    # Depth-first search to check chain connectivity and linearity
    while stack:
        current_atom, previous_atom = stack.pop()
        if current_atom.GetIdx() in visited:
            continue
        visited.add(current_atom.GetIdx())

        # Check for non-single bonds
        for bond in current_atom.GetBonds():
            neighbor = bond.GetOtherAtom(current_atom)
            if neighbor.GetIdx() == (previous_atom.GetIdx() if previous_atom else None):
                continue
            if neighbor.GetIdx() not in visited:
                if neighbor.GetAtomicNum() == 6 and bond.GetBondType() == Chem.BondType.SINGLE:
                    stack.append((neighbor, current_atom))

    # Count total carbon atoms in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Ensure the visited carbon atoms make up the entire chain
    if len(visited) < carbon_count:
        return False, "Carbon chain is branched or has unvisited sections"

    return True, "Molecule is a straight-chain saturated fatty acid"

# Example for testing the implementation
example_smiles = "CCCCCCC(O)=O"  # Heptanoic acid
result, reason = is_straight_chain_saturated_fatty_acid(example_smiles)
print(f"{example_smiles}: {result} ({reason})")