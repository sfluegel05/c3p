"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid has one or more alkyl substituents and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (R-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, "Insufficient carbon count for alkyl chain"

    # Look for branching - allow various forms of branching including methyl branches
    branch_pattern = Chem.MolFromSmarts("[CX4;!R][CX4;!R,CX3]([CX4;!R,CX3,H])[CX4;!R,CX3,H]")
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No branch points found"

    # Check for excessive ring structures, which might not be typical for branched chain fatty acids
    if mol.GetRingInfo().NumRings() > 1:
        return False, "Too many ring structures for a typical branched-chain fatty acid"

    # Verify the branch complexity isn't too high for typical fatty acid structures
    branch_complexity = 0
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:  # Assumes branches will have more than 2 neighbors
            branch_complexity += atom.GetDegree()
    if branch_complexity > 15:
        return False, "Structure too complex for a typical branched-chain fatty acid"

    return True, "Contains a carbon chain with branch points and a carboxylic acid group with acceptable complexity"