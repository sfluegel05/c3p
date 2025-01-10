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
    
    # Define the corrin nucleus SMARTS pattern
    # This pattern should represent a more realistic approximation of the corrin macrocycle
    # Full understanding of pyrrole and its reduced forms, along with connecting carbons is vital
    # Here we attempt to match four pyrrole rings connected in a specific, macrocyclic fashion
    # Due to complexity, this is an educated guess and may need refinement for accuracy in real databases
    corrin_core_pattern = Chem.MolFromSmarts("C1=C(C2=N[C@@H](C)C=C3C=NC(C[C@@H]4C[NH]C(=C1)[C@@H]4C3=CC2)C)C")
    
    # Match molecule to the corrin pattern
    if mol.HasSubstructMatch(corrin_core_pattern):
        return True, "Contains a corrin nucleus pattern"
    
    return False, "Does not contain a corrin nucleus pattern"