"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline or its hydrogenated derivative based on its SMILES string.
    A beta-carboline is defined as any pyridoindole containing a beta-carboline skeleton and their hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline or hydrogenated derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define beta-carboline skeleton SMARTS pattern with variable bond orders
    beta_carboline_smarts = '[nH]1cc2ccccc2c3cccnc13'
    beta_carboline_mol = Chem.MolFromSmarts(beta_carboline_smarts)
    if beta_carboline_mol is None:
        return False, "Error generating beta-carboline query molecule"

    # Perform substructure search
    if not mol.HasSubstructMatch(beta_carboline_mol):
        return False, "Beta-carboline skeleton not found"

    return True, "Contains beta-carboline skeleton or its hydrogenated derivative"