"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine contains an alkyl group attached to an aromatic ring and an amine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # A specific SMARTS pattern for an aralkylamine
    # An aromatic ring bonded to an alkyl chain which then links to an amine
    aralkylamine_pattern = Chem.MolFromSmarts("a[A;R0][CX4][NX3;H2,H1,H0]")

    if not mol.HasSubstructMatch(aralkylamine_pattern):
        return False, "Structure doesn't match aralkylamine pattern"

    return True, "Contains an aromatic ring linked via alkyl to an amine group"