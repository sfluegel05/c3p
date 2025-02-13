"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A BCFA is characterized by the presence of a carboxylic acid group and one or more branches
    in its hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a BCFA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Locate branching - carbons with more than two heavy atom (non-hydrogen) bonds
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            neighbors = atom.GetNeighbors()
            if len([nbr for nbr in neighbors if nbr.GetAtomicNum() > 1]) > 2:
                # Found a branched carbon
                return True, "Contains a branched hydrocarbon chain with carboxylic acid group"

    return False, "No branching carbon chain found"

# Examples:
smiles_examples = [
    "OC(=O)C=C(C)C",  # Example of branched-chain fatty acid
    "CCCCCCCCCCCCCCCCCC(=O)O",  # Example of straight-chain fatty acid
]

results = {sm: is_branched_chain_fatty_acid(sm) for sm in smiles_examples}
for sm, result in results.items():
    print(f"SMILES: {sm} -> {result}")