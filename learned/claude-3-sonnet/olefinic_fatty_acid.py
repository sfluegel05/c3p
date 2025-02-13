"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: CHEBI:51770 olefinic fatty acid
An olefinic fatty acid is any fatty acid containing at least one C=C double bond.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.

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

    # Look for carboxylic acid group (-C(=O)O-) at the end of a linear chain
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O[C;H2]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No valid carboxylic acid group found"

    # Look for carbon-carbon double bonds (C=C) in a linear chain
    double_bond_pattern = Chem.MolFromSmarts("[C;H2]=C[C;H2]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if not double_bond_matches:
        return False, "No carbon-carbon double bonds found in a linear chain"

    # Check for long linear carbon chain (fatty acid)
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() <= 2)
    if n_carbons < 6 or n_carbons > 24:  # Typical fatty acid range
        return False, "Carbon chain length outside the typical range for fatty acids"

    # Additional checks for specific double bond configurations or other criteria

    return True, "Contains a carboxylic acid group and at least one carbon-carbon double bond in a linear chain"