"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone with a germacrane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for germacrane skeleton (10-membered carbon ring)
    germacrane_pattern = Chem.MolFromSmarts('C1CCCCC[CH]CCC1')  # simplified carbon ring structure
    if not mol.HasSubstructMatch(germacrane_pattern):
        return False, "No germacrane skeleton found"
    
    # Look for lactone group (-C(=O)O-)
    lactone_pattern = Chem.MolFromSmarts('C(=O)O')
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if len(lactone_matches) < 1:
        return False, "No lactone group found"
    
    # Check for the presence of typical sesquiterpene traits like double bonds/epoxides
    # Note: This is a simplification; exact stereochemistry consideration requires more complex checks
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if double_bond_count < 2:
        return False, "Insufficient number of double bonds for sesquiterpene lactone"

    # (Optional) Further detailed stereochemical checks and substituents specific to known germacranolides

    return True, "Contains germacrane skeleton with lactone group, typical of germacranolides"