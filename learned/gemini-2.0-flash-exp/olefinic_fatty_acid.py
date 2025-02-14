"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Any fatty acid containing at least one C=C double bond.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid is a fatty acid with at least one C=C double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for at least one C=C double bond
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No C=C double bond found"

    # Check for a reasonable number of carbons to be considered a fatty acid (at least 4)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if num_carbons < 4:
         return False, "Too few carbons to be a fatty acid."

    return True, "Contains a carboxylic acid group and at least one C=C double bond and has 4 or more carbons"