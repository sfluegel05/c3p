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

    # Look for evidence of branching by checking for non-linear shapes
    # Consider a simple search for carbon with more than two bonds as potential branch points
    branch_pattern = Chem.MolFromSmarts("[C]([C])([C])")
    if not mol.HasSubstructMatch(branch_pattern) and not Chem.MolFromSmarts("[CX4](C)C"):
        return False, "No branch points found"

    # Check for ring structures, considering that some BCFAs can indeed have small rings like cyclopropyl
    # We will allow BCFAs to have one or two rings as this is not uncommon in complex natural BCFAs
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 2:  
        return False, "Too many ring structures for a typical branched-chain fatty acid"

    return True, "Contains a carbon chain with branch points and a carboxylic acid group with acceptable complexity"