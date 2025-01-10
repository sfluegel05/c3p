"""
Classifies: CHEBI:133249 saturated fatty aldehyde
"""
from rdkit import Chem

def is_saturated_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a saturated fatty aldehyde based on its SMILES string.
    A saturated fatty aldehyde has an aldehyde group with a saturated (no carbon-carbon double bonds) carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aldehyde group pattern: O=C and ensuring the terminal carbon is present
    aldehyde_pattern = Chem.MolFromSmarts("[CX3](=O)[!$([#1])]")  # Focus on carbon doubly bonded to oxygen
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"
    
    # Ensure no carbon-carbon double bonds: restrict to single bonds only
    unsat_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(unsat_pattern):
        return False, "Contains carbon-carbon unsaturation"

    # Check if it's a fatty aldehyde - usually long carbon chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Too few carbons to be a fatty aldehyde"

    return True, "Contains aldehyde group with a saturated carbon chain"