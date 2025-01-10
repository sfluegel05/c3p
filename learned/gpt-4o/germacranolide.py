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
    
    # More generalized germacrane skeleton pattern (can be refined iteratively)
    germacrane_macrocycle_pattern = Chem.MolFromSmarts('C1C=CCC=CC=CC1')
    if not mol.HasSubstructMatch(germacrane_macrocycle_pattern):
        return False, "No germacrane skeleton found (broad pattern)"
    
    # Look for a lactone moiety (OC1=O in ring context)
    lactone_pattern = Chem.MolFromSmarts('O=C1OC=CC=CC1')
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group in a cyclic context found"

    # Number of double bonds - typical region variability should allow for different locations
    double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE]
    if len(double_bonds) < 2:
        return False, "Insufficient double bonds for sesquiterpene nature"

    # Optional: Further checks on stereochemistry if very necessary, or consider relaxed rules
    # These are often highly variable in natural systems, so use with caution

    return True, "Contains a germacrane skeleton with a macrocyclic lactone, typical of germacranolides"