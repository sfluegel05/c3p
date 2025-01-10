"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    A beta-carboline is defined as any pyridoindole containing a beta-carboline skeleton 
    and their hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an enhanced SMARTS pattern for beta-carboline core (pyridoindole structure)
    # capturing more variants found in example structures
    beta_carboline_patterns = [
        Chem.MolFromSmarts("c12c([nH]c3ccccc13)nccc2"),  # Basic beta-carboline core
        Chem.MolFromSmarts("c12c([nH]c3cccc[nH]13)nccc2"),  # Dihydro-beta-carboline
        Chem.MolFromSmarts("c12c([nH]c3ccccc13)N[CH2]C2"),  # Hydrogenated and substituted
        Chem.MolFromSmarts("c12c3ccccc3[nH]c1[nH]cc2"),  # Fully aromatic with indole
        Chem.MolFromSmarts("c12c3ccccn3c1[nH]c[nH]2"),  # Variably substituted
        Chem.MolFromSmarts("n1c2c(cccc2)c2c1n[nH]c2"),  # Pyridoindole with variations
    ]

    # Check for matching any of the expanded beta-carboline patterns
    for pattern in beta_carboline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-carboline core structure"

    return False, "Does not contain beta-carboline core structure"