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
    
    # SMARTS pattern for detecting a corrin-like structure
    # This pattern tries to capture the four pyrrole rings connectivity and macrocyclic nature
    # Example SMARTS for pyrrole: c1cc[nH]c1; 
    # Modification here is a task-specific guess and must be refined
    corrin_pattern = Chem.MolFromSmarts(
        "[nH]1cc[cH][cH]1-[cH]2[cH][cH][c][cH]2"  # simplified and abstract; more detail would be needed
    )  
    
    # Match the molecule to the corrin pattern
    if mol.HasSubstructMatch(corrin_pattern):
        return True, "Contains a corrin nucleus pattern"

    return False, "Does not contain a corrin nucleus pattern"