"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on a germacrane skeleton.

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
    
    # Proposed germacrane ring + lactone pattern with a broader search
    # This is a simplification, capturing needed structural pieces
    germacrane_macrocycle_pattern = Chem.MolFromSmarts('C1(C2)C=CCCC(=CC1)C(=O)OC2')
    if not mol.HasSubstructMatch(germacrane_macrocycle_pattern):
        return False, "No germacrane macrocycle with lactone structure found"
    
    # Check for presence of at least two double bonds; germacrane typically has multiple double bonds
    double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE]
    if len(double_bonds) < 2:
        return False, "Insufficient double bonds for sesquiterpene nature"

    # Check for stereochemistry features: comment this if it's known to be variable in context
    # Consider specific known patterns for stereocenters if necessary

    return True, "Contains a germacrane skeleton with a macrocyclic lactone, typical of germacranolides"