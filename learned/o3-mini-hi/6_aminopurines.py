"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies chemical entities for the class 6-aminopurines.
Definition: Any compound having 6-aminopurine (adenine) as part of its structure.
This version uses two adenine substructure templates to be robust against different representations.
"""

from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines whether a molecule contains a 6-aminopurine (adenine) substructure.
    To account for different representations, two templates are checked.
    
    Args:
        smiles (str): SMILES string of the compound.
    
    Returns:
        bool: True if the molecule contains an adenine (6-aminopurine) substructure, False otherwise.
        str: A reason explaining the decision.
    """
    # Parse the input SMILES to create an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the adenine (6-aminopurine) substructure templates.
    # Template 1: common representation in many CoA derivatives (using aromatic notation).
    adenine_template1 = Chem.MolFromSmarts("n1cnc(N)c2ncnc12")
    if adenine_template1 is None:
        return False, "Failed to create adenine template 1"
    
    # Template 2: alternative aromatic representation of adenine.
    adenine_template2 = Chem.MolFromSmarts("c1nc2c(n1)nc(nc2)N")
    if adenine_template2 is None:
        return False, "Failed to create adenine template 2"
    
    # Check if the molecule contains either adenine substructure.
    if mol.HasSubstructMatch(adenine_template1) or mol.HasSubstructMatch(adenine_template2):
        return True, "Contains 6-aminopurine (adenine) substructure"
    else:
        return False, "6-aminopurine substructure not found"

# Example usage:
if __name__ == "__main__":
    # Test with Triacanthine, which is known to contain the adenine moiety.
    test_smiles = "CC(C)=CCn1cnc(N)c2ncnc12"
    result, reason = is_6_aminopurines(test_smiles)
    print(f"Result: {result}\nReason: {reason}")