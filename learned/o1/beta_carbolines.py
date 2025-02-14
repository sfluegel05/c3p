"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: Beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline or its hydrogenated derivative based on its SMILES string.
    Beta-carbolines are pyridoindoles containing a beta-carboline skeleton and their hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-carboline core pattern (pyrido[3,4-b]indole)
    # This SMARTS pattern represents the fused tricyclic ring system
    beta_carboline_smarts = """
    [$([nH]),n]1ccc2c1[nH,c]c3cccc[n,c]23
    """
    beta_carboline_pattern = Chem.MolFromSmarts(beta_carboline_smarts)
    if beta_carboline_pattern is None:
        return False, "Failed to create beta-carboline pattern"

    # Check for substructure match
    if mol.HasSubstructMatch(beta_carboline_pattern):
        return True, "Contains beta-carboline core"
    else:
        return False, "Does not contain beta-carboline core"