"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester is any carbonate in which the hydrogens of carbonic acid
    have been replaced by organyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for carbonate ester functional group
    # Carbon atom double-bonded to oxygen and single-bonded to two oxygens,
    # each oxygen single-bonded to a non-hydrogen atom (organyl group)
    carbonate_pattern = Chem.MolFromSmarts("[CX3](=O)(O[!H])(O[!H])")
    if carbonate_pattern is None:
        return False, "Invalid SMARTS pattern for carbonate ester"
    
    # Search for the carbonate ester functional group
    matches = mol.GetSubstructMatches(carbonate_pattern)
    if matches:
        return True, "Contains carbonate ester functional group"
    else:
        return False, "Does not contain carbonate ester functional group"
    
    
__metadata__ = {
    'chemical_class': {
        'name': 'carbonate ester',
        'definition': 'Any carbonate that is carbonic acid in which the hydrogens have been replaced by organyl groups.',
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
    },
}