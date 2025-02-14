"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: CHEBI:51714 enone
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha,beta-unsaturated ketone of general formula
    R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) ≠ H) in which the C=O function is
    conjugated to a C=C double bond at the alpha,beta position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for enone
    # [C;H0]=[C]-[C](=O)-[C;!H0]
    # This matches a double bond between two carbons (alpha and beta),
    # followed by a single bond to a carbonyl group (ketone),
    # where the carbonyl carbon has no attached hydrogens (not an aldehyde)
    enone_smarts = '[C;H0]-[C]=[C]-C(=O)[C;H0]'

    enone_pattern = Chem.MolFromSmarts(enone_smarts)
    if enone_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for enone pattern in molecule
    matches = mol.GetSubstructMatches(enone_pattern)
    if matches:
        return True, "Contains enone group (alpha,beta-unsaturated ketone)"
    else:
        return False, "Does not contain enone group"