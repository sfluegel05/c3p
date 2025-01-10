"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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

    # Check for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Check for unsaturation (double/triple bonds) or branches
    query = rdqueries.AtomNumEqualsQueryAtom(6) # carbon atoms
    unsaturations = Chem.MolFromSmarts("[C]=[C]")  # unsaturated bonds
    if mol.HasSubstructMatch(unsaturations):
        return False, "Contains unsaturated (double/triple) bonds"

    # Check for branches - all carbon atoms must be in a single line (acyclic unbranched carbon chain)
    # We iterate and check carbon neighbors
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6: # Carbon
            # If carbon is connected to more than 2 other carbons, it indicates branching
            carbon_neighbors = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            if carbon_neighbors > 2:
                return False, "Contains branched carbon chain"

    return True, "Contains a linear, saturated carbon chain with a terminal carboxylic acid group"