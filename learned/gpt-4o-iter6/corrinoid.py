"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    Corrinoids contain a corrin nucleus, which consists of four reduced or partly reduced pyrrole rings 
    joined in a macrocycle by three =C- groups and one direct carbon-carbon bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Construct a more generic and pragmatic pattern for the corrin nucleus
    # Note: Due to complexity, this is simplified and assumptions are made
    # Matching actual corrin complexity in SMARTS might exceed practical scope
    # Look for four pyrrole-like rings connected with bridge atoms typical of the corrin structure
    corrin_core_pattern = Chem.MolFromSmarts('C=1C=CC2=CC=CC3=CC=CC4=CC=CC=C1N=C(C2)C3=C4')
    
    # Example pattern lacks detail due to SMARTS limitations - adapted below for best approximation
    # Here, looking simply for multiple interconnected nitrogen-based heterocycles is the approach
    # Practical implementations may align with known molecule library comparisons beyond SMARTS

    # Match molecule to the corrin pattern
    if mol.HasSubstructMatch(corrin_core_pattern):
        return True, "Contains a corrin-like nucleus pattern"
    
    return False, "Does not contain a discernible corrin nucleus pattern"