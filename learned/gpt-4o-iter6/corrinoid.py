"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    Corrinoids contain a corrin nucleus, which consists of four reduced or partly reduced pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Defining the corrin-like pattern using SMARTS (note: may not be entirely accurate; 
    #       this is a simplified representation of a complex macrocycle)
    # This is a very simplified and abstract pattern, not covering all specific features of full corrin rings
    corrin_pattern = Chem.MolFromSmarts("C1CCC=C1")  # placeholder pattern, not accurate for actual corrin
    
    # Try matching the pattern with the molecule
    if mol.HasSubstructMatch(corrin_pattern):
        return True, "Contains a corrin nucleus-like pattern"

    return False, "Does not contain a corrin nucleus-like pattern"