"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    if len(carboxylic_matches) != 1:
        return False, "No carboxylic acid group found or multiple groups present"

    # Ensure the terminal carbon is part of the carboxylic acid
    terminal_carbon = carboxylic_matches[0][0]
    # Ensure it's the end of the chain
    neighbors = mol.GetAtomWithIdx(terminal_carbon).GetNeighbors()
    if any(neighbor.GetAtomicNum() == 6 for neighbor in neighbors if neighbor.GetIdx() != terminal_carbon + 1):
        return False, "Carboxylic acid group not at the end of a chain"

    # Ensure no double or triple bonds in carbon chain
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False, "Presence of non-single C-C bonds indicates unsaturation"

    # Check for hydroxy groups, allow at most one (apart from carboxylic acid O)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]C")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)

    allowed_hydroxy_count = 1
    if len(hydroxy_matches) > allowed_hydroxy_count:
        return False, f"Too many hydroxy groups: found {len(hydroxy_matches)}, allowed at most {allowed_hydroxy_count}"

    # Ensure it's a straight chain
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]

    # Create a linear chain from terminal carbon and check connectivity
    current = terminal_carbon
    chain_atoms = set()
    for _ in range(len(carbon_atoms)):
        chain_atoms.add(current)
        neighbors = [
            neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(current).GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in chain_atoms
        ]
        if len(neighbors) != 1:
            return False, "Branches or side chains detected; not a straight chain"
        current = neighbors[0]

    if len(chain_atoms) != len(carbon_atoms):
        return False, "Not a continuous straight chain; branches or side chains detected"

    return True, "Molecule is a straight-chain saturated fatty acid"