"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid is an unsaturated fatty acid with at least one double bond and a terminal carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for the presence of a carbon-carbon double bond
    olefinic_double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(olefinic_double_bond_pattern):
        return False, "No olefinic double bonds found"
    
    # Check that there are at least 8 carbons to consider it a fatty acid
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 8:
        return False, "Too few carbons for a fatty acid"
    
    # If passed all criteria, the molecule is classified as an olefinic fatty acid
    return True, "Contains olefinic double bonds and a carboxylic acid group"

# Sample usage of the function:
#print(is_olefinic_fatty_acid("O(O)C(CCCCC)/C=C/C=C\C/C=C\CCCCC(O)=O"))  # Example SMILES