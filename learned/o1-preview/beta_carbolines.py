"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline or its hydrogenated derivative based on its SMILES string.
    A beta-carboline is a pyridoindole containing a beta-carboline skeleton and their hydrogenated derivatives.

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

    # Define the beta-carboline core skeleton using SMARTS
    # This pattern represents the fused pyridine and indole ring system
    beta_carboline_smarts = """
    [#7H]1c2ccccc2c2nccccc12   # beta-carboline skeleton
    """
    # Remove comments and whitespace in SMARTS
    beta_carboline_smarts = beta_carboline_smarts.strip().split('#')[0].strip()

    beta_carboline_pattern = Chem.MolFromSmarts(beta_carboline_smarts)
    if beta_carboline_pattern is None:
        return False, "Invalid SMARTS pattern for beta-carboline skeleton"

    # Check for beta-carboline skeleton
    if mol.HasSubstructMatch(beta_carboline_pattern):
        return True, "Contains beta-carboline skeleton"
    else:
        # Allow for hydrogenated derivatives by allowing non-aromatic atoms and bonds
        hydrogenated_smarts = """
        [#7H]1[c;R][c;R][c;R][c;R][c;R]2[c;R]1[c;R][n;R][c;R][c;R][c;R]2    # hydrogenated beta-carboline skeleton
        """
        hydrogenated_smarts = hydrogenated_smarts.strip().split('#')[0].strip()
        hydrogenated_pattern = Chem.MolFromSmarts(hydrogenated_smarts)
        if hydrogenated_pattern is None:
            return False, "Invalid SMARTS pattern for hydrogenated beta-carboline"

        if mol.HasSubstructMatch(hydrogenated_pattern):
            return True, "Contains hydrogenated beta-carboline skeleton"
        else:
            return False, "Does not contain beta-carboline skeleton"