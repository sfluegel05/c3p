"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies chemical entities for the class 6-aminopurines.
Definition: Any compound having 6-aminopurine (adenine) as part of its structure.
This improved version not only searches for a SMARTS match but also checks that the match,
when extracted from the molecule, is chemically equivalent to a canonical adenine fragment.
"""

from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) substructure.
    This version uses an exact fragment matching strategy:
      1. An adenine template (6-aminopurine) is defined from SMILES.
      2. The molecule is searched for a substructure match.
      3. For every match, we extract that fragment and generate its canonical SMILES.
         If it exactly equals the canonical SMILES of adenine, then we classify the input
         as containing an unmodified adenine unit.
    
    Args:
        smiles (str): The SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule contains an adenine fragment, False otherwise.
        str: A reason explaining the classification.
    """
    # Parse the input SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the adenine (6-aminopurine) template and compute its canonical SMILES.
    adenine_template = Chem.MolFromSmiles("n1cnc2c(N)ncnc12")
    if adenine_template is None:
        return False, "Failed to create adenine template"
    canonical_adenine = Chem.MolToSmiles(adenine_template, canonical=True)
    
    # Define the SMARTS pattern using the same pattern as used for the template.
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if adenine_pattern is None:
        return False, "Failed to define adenine SMARTS pattern"
    
    # Get substructure matches for adenine in the molecule.
    matches = mol.GetSubstructMatches(adenine_pattern)
    if not matches:
        return False, "6-aminopurine substructure not found"
    
    # For every match, extract the fragment and compare its canonical SMILES to adenine.
    # We use Chem.PathToSubmol to extract the submolecule defined by the match.
    for match in matches:
        fragment = Chem.PathToSubmol(mol, match)
        fragment_smiles = Chem.MolToSmiles(fragment, canonical=True)
        if fragment_smiles == canonical_adenine:
            return True, "Contains 6-aminopurine (adenine) substructure"
    
    # If none of the matches exactly equal adenine, we assume the structure has modifications
    # that prevent a strict assignment to the 6-aminopurine class.
    return False, "Found a purine-like substructure but not an unmodified 6-aminopurine"
    
# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: Triacanthine.
    test_smiles = "CC(C)=CCn1cnc(N)c2ncnc12"
    result, reason = is_6_aminopurines(test_smiles)
    print(f"Result: {result}\nReason: {reason}")