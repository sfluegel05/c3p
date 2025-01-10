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

    # Look for branch points; carbon atoms with more than one carbon substituent
    branch_pattern = Chem.MolFromSmarts("[CX4]([CX4])[CX4]")
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No alkyl branching found"
    
    # Allow for small rings if present; ensure they are limited in number
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 2:
        return False, "Too many ring structures for a typical branched-chain fatty acid"

    # Check that branches are not too complex (exclude non-alkyl groups)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [7, 8, 16]:  # Exclude nitrogen, oxygen, sulfur in branches
            if any(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors()):
                return False, "Presence of non-alkyl groups in branches"

    return True, "Contains a carbon chain with branch points and a carboxylic acid group with acceptable complexity"