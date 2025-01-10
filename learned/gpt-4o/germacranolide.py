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
    
    # Detailed pattern matching for the germacrane 10-membered carbon ring.
    germacrane_pattern = Chem.MolFromSmarts('C1C=CCC2OC(=O)C=CC(=O)C1C2') 
    if not mol.HasSubstructMatch(germacrane_pattern):
        return False, "No specific germacrane skeleton found"

    # Look for lactone group specifically within macrocyclic context
    lactone_pattern = Chem.MolFromSmarts('C1OC(=O)C1')
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if len(lactone_matches) < 1:
        return False, "No lactone group in macrocyclic ring found"

    # Count double bonds to enforce typical sesquiterpene traits
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if double_bond_count < 2:
        return False, "Insufficient number of double bonds for sesquiterpene lactone"

    # Check stereochemistry if needed, this can be expanded with more specific patterns for known stereochemistry
    # Consider known substituents common in germacranolides

    return True, "Contains structured germacrane skeleton with lactone group, typical of germacranolides"