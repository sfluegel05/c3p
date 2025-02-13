"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine in which the alkyl group is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an aralkylamine, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for an aromatic group
    # Ensure the molecule contains at least one aromatic ring
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("a")):
        return False, "No aromatic group found"

    # Ensure the molecule contains an amine group (primary, secondary, or tertiary, but not quaternary)
    amine_pattern = Chem.MolFromSmarts("[NX3;!$([N+])]")  # Amine without quaternary nitrogen
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"

    # Check for aralkyl-amine connections
    # SMARTS pattern for an aromatic group connected via an alkyl chain to an amine
    aralkylamine_pattern = Chem.MolFromSmarts("[a]C[NX3;!$([N+])]")  # Aromatic-C-alkyl-N
    if not mol.HasSubstructMatch(aralkylamine_pattern):
        return False, "Aromatic group not properly attached through an alkyl to amine"

    return True, "Contains an alkylamine with aromatic substitution"

# Testing the function with examples should follow to ensure it meets the criteria correctly.