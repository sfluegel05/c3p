"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.

    A methyl-branched fatty acid is defined as any branched-chain fatty acid
    containing methyl branches only and having a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for methyl branches on acyl chain
    # Adjust the pattern to match more general methyl branching (smallest alkyl branch possible)
    methyl_branch_pattern = Chem.MolFromSmarts("[C;!R][CH2][C](C)O")
    if not mol.HasSubstructMatch(methyl_branch_pattern):
        return False, "Molecule does not contain methyl branches on a fatty acid backbone"
    
    # Estimate chain length by counting carbons and ensure it's within expected range for fatty acids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Carbon chain too short for traditional fatty acid"
    
    return True, "Molecule is a methyl-branched fatty acid"