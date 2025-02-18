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
    # Core structure with aromatic rings
    core_pattern = Chem.MolFromSmarts("[nH]1c2ccccc2c2cnccc12")
    if core_pattern is not None and mol.HasSubstructMatch(core_pattern):
        return True, "Contains aromatic beta-carboline core structure"

    # Check for hydrogenated derivatives (single bonds in the rings)
    # Flexible pattern allowing single bonds
    flexible_pattern = Chem.MolFromSmarts("[nH]1-C-C-C2-C1-C1-C-C-C3-C1-N-C-C-C2-C3")
    if flexible_pattern is not None and mol.HasSubstructMatch(flexible_pattern):
        return True, "Contains hydrogenated beta-carboline core"

    # If no matches, return False
    return False, "No beta-carboline core structure detected"