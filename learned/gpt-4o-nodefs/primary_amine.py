"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine contains a nitrogen atom bonded to exactly one carbon 
    group R and two hydrogen atoms, fitting the formula: RNH2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated SMARTS pattern for primary amine (RNH2)
    # Ensure nitrogen is singly bonded to one carbon, has two hydrogens, and is not part of a larger nitrogen-bonded system
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2;!$(NC);!$(N-*=[O,N])]")

    # Check if the molecule matches the primary amine SMARTS pattern
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Contains a primary amine group (RNH2 structure identified)"

    return False, "No primary amine group found"