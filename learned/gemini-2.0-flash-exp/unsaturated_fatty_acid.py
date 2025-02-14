"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid contains at least one C=C or C#C bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an unsaturated fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for C=C and C#C bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")

    # Check for the presence of either C=C or C#C bonds
    if mol.HasSubstructMatch(double_bond_pattern) or mol.HasSubstructMatch(triple_bond_pattern):
        return True, "Contains at least one C=C or C#C bond"
    else:
         return False, "Does not contain any C=C or C#C bond"