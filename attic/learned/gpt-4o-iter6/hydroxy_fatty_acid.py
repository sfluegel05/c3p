"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is a fatty acid with one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find the carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Find hydroxy groups, ensuring they are not aromatic (i.e., on sp2 carbon)
    hydroxy_pattern = Chem.MolFromSmarts("[CX4;!$(C=[O,N])][OH]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 1:
        return False, "No hydroxy groups found"

    # Count number of non-ring carbons to identify linear chain
    linear_carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing())

    # Allow shorter chains (minimum 4) recognizing short hydroxy fatty acids
    if linear_carbon_count < 4:
        return False, "Carbon chain too short for a typical fatty acid"

    return True, "Contains one or more hydroxy groups and a carboxylic acid group with sufficient carbon chain length"