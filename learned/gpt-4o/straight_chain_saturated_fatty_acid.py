"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    A straight-chain saturated fatty acid has a linear carbon chain ending in a carboxylic acid group
    without any side chains, but may contain hydroxy groups at linear positions.

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

    # Pattern for linear saturated carbon chain with terminal carboxylic acid group
    fatty_acid_pattern = Chem.MolFromSmarts("C(CCCCCCCCCCCCCCCCC(O)=O)CCC")

    # Check for the pattern within the molecule
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "Does not match pattern for a straight-chain saturated fatty acid"

    # Scan atoms to ensure no branching, calculate carbon chain length
    carbon_chain_length = 0
    last_carbon_idx = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            carbon_chain_length += 1
            neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
            # Check that carbons bonded to no more than 2 other carbons
            if neighbors.count('C') > 2:
                return False, f"Contains branching at carbon {atom.GetIdx()}"
            # Also fetch last carbon index
            last_carbon_idx = atom.GetIdx()

    # Verify terminal position ends with carboxylic acid and hydroxyl can exist
    term_carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(term_carboxyl_pattern):
        return False, "Does not end in a carboxylic acid group"

    # Allow hydroxyls internally (without branching)
    hydroxyl_groups = Chem.MolFromSmarts("[CH3,CH2]CO")
    if not all(mol.HasSubstructMatch(hydroxyl_groups) or neighbor.GetSymbol() == 'O' for neighbor in mol.GetAtoms()):
        return False, "Contains disallowed functional groups"

    return True, f"Matches straight-chain saturated fatty acid with length {carbon_chain_length}"