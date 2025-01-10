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
    
    # Check for terminal aldehyde group pattern: R-CHO
    terminal_aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CH2]")
    if not mol.HasSubstructMatch(terminal_aldehyde_pattern):
        return False, "Terminal aldehyde group not found"

    # Ensure only one aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if len(mol.GetSubstructMatches(aldehyde_pattern)) != 1:
        return False, "More than one aldehyde group present or not at terminal position"

    # Ensure no carbon-carbon double bonds
    unsaturation_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(unsaturation_pattern):
        return False, "Contains carbon-carbon unsaturation"

    # Check for a minimum of 5 consecutive carbons, including aldehyde carbon
    carbon_chain_pattern = Chem.MolFromSmarts("[CH2]" * 4 + "[CX3H1](=O)[CH2]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Insufficient chain length including aldehyde carbon (minimum 5 carbons)"

    # Ensure no other oxygen-containing functional groups present
    other_oxygen_pattern = Chem.MolFromSmarts("[OX2H1,OX2H0]")
    if mol.HasSubstructMatch(other_oxygen_pattern):
        return False, "Additional oxygen-containing functional groups detected"

    # Ensure no rings - fatty aldehydes are non-cyclic
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Contains rings, not a linear/branched carbon chain"

    return True, "Contains terminal aldehyde group with a saturated carbon chain"