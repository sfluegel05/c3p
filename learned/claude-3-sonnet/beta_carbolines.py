"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies: CHEBI:33243 beta-carbolines

A beta-carboline is defined as any pyridoindole containing a beta-carboline skeleton
and their hydrogenated derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_carboline(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    except ValueError:
        return False, "Invalid SMILES string"
    
    # Check for pyridoindole core
    pyridoindole_pattern = Chem.MolFromSmarts("c1nc2ccccc2c3cccnc13")
    if not mol.HasSubstructMatch(pyridoindole_pattern):
        return False, "Does not contain pyridoindole core"
    
    # Check for beta-carboline skeleton
    beta_carboline_pattern = Chem.MolFromSmarts("c1nc2c3ccccc3nc2c4cccnc14")
    beta_carboline_matches = mol.GetSubstructMatches(beta_carboline_pattern, useChirality=False)
    
    # Check for hydrogenated derivatives
    hydrogenated_pattern = Chem.MolFromSmarts("c1nc2c3ccccc3nc2c4ccccn14")
    hydrogenated_matches = mol.GetSubstructMatches(hydrogenated_pattern, useChirality=False)
    
    if beta_carboline_matches or hydrogenated_matches:
        return True, "Contains beta-carboline skeleton or hydrogenated derivative"
    else:
        return False, "Does not contain beta-carboline skeleton or hydrogenated derivative"