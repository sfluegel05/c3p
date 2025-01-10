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

    # Define the beta-carboline core structure using SMILES (pyrido[3,4-b]indole)
    beta_carboline_smiles = 'c1cc2[nH]c3ccccc3c2cn1'  # Beta-carboline core
    beta_carboline_mol = Chem.MolFromSmiles(beta_carboline_smiles)

    if beta_carboline_mol is None:
        return False, "Invalid SMILES for beta-carboline core structure"

    # Check for beta-carboline skeleton
    if mol.HasSubstructMatch(beta_carboline_mol):
        return True, "Contains beta-carboline skeleton"
    else:
        # Check for hydrogenated derivatives
        # Generalized pattern allowing for saturation in the rings
        hydrogenated_smarts = '[nH]1c2cccc(c2c2ccccn12)'  # Generalized pattern
        hydrogenated_pattern = Chem.MolFromSmarts(hydrogenated_smarts)

        if hydrogenated_pattern is None:
            return False, "Invalid SMARTS for hydrogenated beta-carboline core"

        if mol.HasSubstructMatch(hydrogenated_pattern):
            return True, "Contains hydrogenated beta-carboline skeleton"
        else:
            return False, "Does not contain beta-carboline skeleton"