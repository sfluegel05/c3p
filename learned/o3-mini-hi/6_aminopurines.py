"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies chemical entities for the class 6-aminopurines.
Definition: Any compound having 6-aminopurine (adenine) as part of its structure.
This version uses an exact fragment matching strategy:
  1. An adenine template (6-aminopurine) is defined from SMILES.
  2. The molecule is searched for a substructure match.
  3. For every match, we extract that fragment and generate its canonical SMILES.
     If it exactly equals the canonical SMILES of adenine, then we classify
     the input as containing an unmodified adenine unit.
"""

from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains an unmodified 6-aminopurine (adenine) substructure.
    The strategy used is:
      1. Create an adenine template from the canonical SMILES "Nc1ncnc2ncnc12".
      2. Search the input molecule for substructure matches to a SMARTS pattern that
         corresponds to adenine.
      3. Extract the matched fragment(s) and generate their canonical SMILES.
         If any exactly match that of adenine, then the molecule is classified as containing
         a 6-aminopurine unit.
         
    Args:
        smiles (str): The SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule contains an unmodified adenine fragment, False otherwise.
        str: A reason explaining the classification.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the adenine (6-aminopurine) template using a canonical SMILES representation.
    adenine_template = Chem.MolFromSmiles("Nc1ncnc2ncnc12")
    if adenine_template is None:
        return False, "Failed to create adenine template"
    
    # Generate canonical SMILES for the adenine template.
    canonical_adenine = Chem.MolToSmiles(adenine_template, canonical=True)
    
    # Define a SMARTS pattern for adenine using the same string.
    adenine_pattern = Chem.MolFromSmarts("Nc1ncnc2ncnc12")
    if adenine_pattern is None:
        return False, "Failed to define adenine SMARTS pattern"
    
    # Get substructure matches for adenine in the input molecule.
    matches = mol.GetSubstructMatches(adenine_pattern)
    if not matches:
        return False, "6-aminopurine substructure not found"
    
    # For every match, extract the fragment and check if it is exactly adenine.
    for match in matches:
        # Chem.PathToSubmol extracts the underlying submolecule corresponding to the match.
        fragment = Chem.PathToSubmol(mol, match)
        fragment_smiles = Chem.MolToSmiles(fragment, canonical=True)
        if fragment_smiles == canonical_adenine:
            return True, "Contains 6-aminopurine (adenine) substructure"
    
    # If no fragment exactly equals the adenine template, then the adenine-like part is modified.
    return False, "Found a purine-like substructure but not an unmodified 6-aminopurine"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: Triacanthine.
    test_smiles = "CC(C)=CCn1cnc(N)c2ncnc12"
    result, reason = is_6_aminopurines(test_smiles)
    print(f"Result: {result}\nReason: {reason}")