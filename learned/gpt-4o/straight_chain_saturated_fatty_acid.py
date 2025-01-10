"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    A straight-chain saturated fatty acid has a linear carbon chain ending in a carboxylic acid group,
    without any side chains or unsaturated bonds.

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

    # Define the terminal carboxylic acid group
    carboxylic_acid_terminal_pattern = Chem.MolFromSmarts("C(=O)O")  # A potential improvement by being clear

    # Check if exactly one carboxylic acid group is present and is terminal
    terminal_cidx = None
    idxs = mol.GetSubstructMatches(carboxylic_acid_terminal_pattern)
    if len(idxs) == 1:
        # Get the terminal carbon (attached to the carboxyl group)
        terminal_cidx = idxs[0][0]  # Presume terminal group should be the last to close linearity
    else:
        return False, f"Expected 1 terminal carboxylic acid group, found: {len(idxs)}"
    
    # Check to ensure no branches or unsaturated bonds
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            # Saturation check: must not have double/triple bonds
            for bond in atom.GetBonds():
                if bond.GetBondType().name not in ["SINGLE"]:
                    return False, "Contains unsaturated (non-single) bonds"
            # Branching check: Carbon should only link to two neighbors aside for head/tail carbons
            carbon_neighbors = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            if atom.GetIdx() != terminal_cidx and carbon_neighbors > 2:
                return False, "Contains branched carbon chain"

    # Verify no inappropriate functional groups
    # Ensure only extra functional group allowed is OH when adjacent to the final CH
    allowable_oxygen_subgroup = Chem.MolFromSmarts("[CH]O")  # OH group directly linked to terminal
    if not all(atom.GetSymbol() == "C" or mol.HasSubstructMatch(allowable_oxygen_subgroup) for atom in mol.GetAtoms()):
        return False, "Contains disallowed functional groups, beyond terminal COOH or possible OH"

    # Count additional properties if necessary; e.g., carbon count confirming expected length
    return True, "Contains a linear, saturated carbon chain with a terminal carboxylic acid group"