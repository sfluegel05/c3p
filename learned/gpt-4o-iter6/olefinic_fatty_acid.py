"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid.
    An olefinic fatty acid is characterized by having a carboxylic acid group
    and at least one C=C (carbon-carbon double bond) in an aliphatic chain,
    with certain carbon chain characteristics typical of fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group - allow for common variations
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group or valid variant found"

    # Look for C=C double bonds in an aliphatic environment
    cc_double_bond_pattern = Chem.MolFromSmarts("[C;!R]=[C;!R]")  # Non-ring double bonds
    if not mol.HasSubstructMatch(cc_double_bond_pattern):
        return False, "No non-ring carbon-carbon double bond found"
    
    # Count total number of carbon atoms to verify typical fatty acid length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check for large enough hydrocarbon chain
    if carbon_count < 10:  # Assumes a lower limit for fatty acid
        return False, "Carbon count too low for typical fatty acids"
    
    # Additional checks to exclude non-typical structures:
    # - Avoid compounds with excessive functionalization/perhaps limit oxygen, nitrogen
    # - Exclude known non-fatty structures: heterocycles, etc.
    heteroatom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in (7, 8, 9, 15, 16))
    if heteroatom_count > 5:
        return False, "Excessive heteroatoms for a typical fatty acid"

    return True, "Contains both a carboxylic acid group and at least one aliphatic C=C double bond typical of olefinic fatty acids"