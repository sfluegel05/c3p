"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies chemical entities for the class 6-aminopurines.
Definition: Any compound having 6-aminopurine (adenine) as part of its structure.
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) substructure.
    Adenine is recognized as a purine with an exocyclic amino group on the 6-position.
    Many examples include the SMILES fragment 'n1cnc2c(N)ncnc12'.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if the 6-aminopurine substructure is found, False otherwise.
        str: A reason explaining the classification.
    """
    
    # Parse the input SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for 6-aminopurine (adenine).
    # This pattern captures the purine bicyclic ring with an amino group at position 6.
    adenine_smarts = "n1cnc2c(N)ncnc12"
    adenine_pattern = Chem.MolFromSmarts(adenine_smarts)
    if adenine_pattern is None:
        return False, "Failed to define adenine SMARTS pattern"
    
    # Check if the molecule contains the adenine substructure
    if mol.HasSubstructMatch(adenine_pattern):
        return True, "Contains 6-aminopurine (adenine) substructure"
    else:
        return False, "6-aminopurine substructure not found"
        
# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: Triacanthine.
    test_smiles = "CC(C)=CCn1cnc(N)c2ncnc12"
    result, reason = is_6_aminopurines(test_smiles)
    print(f"Result: {result}\nReason: {reason}")