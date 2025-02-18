"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: beta-carbolines (CHEBI:39176)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
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

    # Define beta-carboline core pattern (pyrido[3,4-b]indole skeleton)
    # This SMARTS accounts for both aromatic and hydrogenated forms
    pattern = Chem.MolFromSmarts("[nH]1c2c(ccc3c1nccc3)cc2")
    if mol.HasSubstructMatch(pattern):
        return True, "Contains beta-carboline core structure"

    # Check for hydrogenated derivatives by allowing single bonds in the rings
    # Alternative pattern with more flexibility in bond types
    flexible_pattern = Chem.MolFromSmarts("[nH]1-*,=c-*,=c-*,=c3c1-*,=n-*,=c-*,=c-3")
    if mol.HasSubstructMatch(flexible_pattern):
        return True, "Contains hydrogenated beta-carboline core"

    # If no matches, return False
    return False, "No beta-carboline core structure detected"