"""
Classifies: CHEBI:133249 saturated fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a saturated fatty aldehyde based on its SMILES string.
    A saturated fatty aldehyde is characterized by an aldehyde group at the end of a saturated (non-cyclic, no carbon-carbon double bonds) carbon chain.

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
    
    # Check for aldehyde group pattern at the terminal position: R-CHO
    aldehyde_pattern = Chem.MolFromSmarts("[#6][CX3H](=O)")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "Aldehyde group not at terminal position"

    # Ensure no carbon-carbon double bonds
    unsaturation_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(unsaturation_pattern):
        return False, "Contains carbon-carbon unsaturation"

    # Check for a minimum of 6 consecutive carbons excluding the aldehyde carbon
    carbon_chain_pattern = Chem.MolFromSmarts("C" * 6)
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Insufficient carbon chain length for a fatty aldehyde"

    # Ensure no rings - fatty alcohols are non-cyclic
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Contains rings, not a linear/branched carbon chain"

    return True, "Contains terminal aldehyde group with a saturated carbon chain"