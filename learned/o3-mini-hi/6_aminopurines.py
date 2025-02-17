"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies chemical entities for the class 6-aminopurines.
Definition: Any compound having 6-aminopurine (adenine) as part of its structure.
This version uses a substructure-search strategy that does not require an exact isolated fragment.
It simply checks whether the adenine (6-aminopurine) core is present somewhere in the molecule.
"""

from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines whether a molecule contains a 6-aminopurine (adenine) substructure.
    Instead of extracting an exact adenine fragment, this function simply checks for
    the presence of the adenine substructure. We use an alternative—but valid—SMILES
    representation for adenine, so that the template is built properly.

    Args:
        smiles (str): SMILES string of the compound

    Returns:
        bool: True if the molecule contains an adenine (6-aminopurine) substructure, False otherwise.
        str: A reason explaining the decision.
    """
    # Parse the input SMILES string to create an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the adenine (6-aminopurine) substructure.
    # Note: Depending on RDKit versions, multiple SMILES can represent adenine.
    # Here we use an alternative representation "NC1=NC2=C(N1)N=CN2".
    adenine_template = Chem.MolFromSmiles("NC1=NC2=C(N1)N=CN2")
    if adenine_template is None:
        return False, "Failed to create adenine template"

    # Instead of trying to extract an exact fragment match (which might fail if the adenine 
    # is fused or has additional substituents), we simply check if the molecule contains
    # the adenine substructure.
    if mol.HasSubstructMatch(adenine_template):
        return True, "Contains 6-aminopurine (adenine) substructure"
    else:
        return False, "6-aminopurine substructure not found"

# Example usage:
if __name__ == "__main__":
    # You can test with one of the provided examples. Here is Triacanthine.
    test_smiles = "CC(C)=CCn1cnc(N)c2ncnc12"
    result, reason = is_6_aminopurines(test_smiles)
    print(f"Result: {result}\nReason: {reason}")