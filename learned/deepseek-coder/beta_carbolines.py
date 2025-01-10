"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: CHEBI:38830 beta-carboline
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    A beta-carboline is a pyrido[3,4-b]indole skeleton, which can be hydrogenated or substituted.

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

    # Define a more flexible beta-carboline core structure pattern
    # This pattern accounts for substitutions and hydrogenations
    beta_carboline_pattern = Chem.MolFromSmarts("[nX2]1c2ccccc2c3c1cc[nX2,cX3]3")
    
    # Check if the molecule contains the beta-carboline core
    if mol.HasSubstructMatch(beta_carboline_pattern):
        return True, "Contains the beta-carboline core structure (pyrido[3,4-b]indole)"
    
    # Define a pattern for hydrogenated derivatives
    # This pattern accounts for fully or partially hydrogenated rings
    hydrogenated_pattern = Chem.MolFromSmarts("[nX2]1c2cccc[cX3,cX2]2c3c1cc[nX2,cX3]3")
    if mol.HasSubstructMatch(hydrogenated_pattern):
        return True, "Contains a hydrogenated beta-carboline core structure"
    
    # Define a pattern for more complex structures with additional rings
    complex_pattern = Chem.MolFromSmarts("[nX2]1c2ccccc2c3c1cc[nX2,cX3]3~*")
    if mol.HasSubstructMatch(complex_pattern):
        return True, "Contains a complex beta-carboline core structure with additional rings or substitutions"
    
    # If no match is found, return False
    return False, "No beta-carboline core structure found"