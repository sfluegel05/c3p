"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid with exactly three double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid or ester group or relevant end group
    carboxylic_acid_or_ester_pattern = Chem.MolFromSmarts("C(=O)[OX1H0]")  # includes ester oxygen or OH
    if not mol.HasSubstructMatch(carboxylic_acid_or_ester_pattern):
        return False, "No carboxylic acid, ester, or similar group found; not a typical trienoic fatty acid structure"

    # Count conjugated C=C double bonds specifically
    # Testing for sequential occurrences common in polyunsaturated conditions
    double_bond_pattern = Chem.MolFromSmarts("C(C=C)C=C(C=C)")  # Seemlessly matches multiple double bonds in sequence
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))

    if double_bond_matches == 1:  # Identifies three double bonds forming a linear chain
        return True, "Contains a typical three double-bond trienoic structure; classified as a trienoic fatty acid"
    else:
        return False, f"Found {double_bond_matches*2} sequence(s) of double bonds; not consistent with classified trienoic fatty acid"

# Note: Adjust the functional group patterns as necessary for specific methyl-branched situations or complex derivations.