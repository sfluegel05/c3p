"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies chemical entities for the class 6-aminopurines.
Definition: Any compound having 6-aminopurine (adenine) as part of its structure.
This version uses the canonical adenine fragment as a substructure query and then 
checks that the extracted fragment exactly matches the adenine skeleton.
"""

from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines whether a molecule contains an unmodified 6-aminopurine (adenine) substructure.
    This function creates an RDKit molecule for adenine (using its canonical SMILES)
    and then performs a substructure search on the input molecule. For every match,
    the corresponding fragment is extracted and its canonical SMILES is compared to that of
    adenine. This helps rule out spurious matches (false positives) that do not exactly contain
    the unmodified adenine core.
    
    Args:
        smiles (str): SMILES string of the compound.
    
    Returns:
        bool: True if the molecule contains an unmodified adenine (6-aminopurine) substructure, False otherwise.
        str: A reason explaining the decision.
    """
    # Parse the input SMILES to create an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Create the adenine molecule using its canonical SMILES.
    adenine_mol = Chem.MolFromSmiles("c1nc2c(n1)nc(nc2)N")
    if adenine_mol is None:
        return False, "Failed to create adenine reference molecule"
    
    # Get the canonical SMILES for adenine (this will be our reference fragment)
    adenine_smiles = Chem.MolToSmiles(adenine_mol, canonical=True)
    
    # Find all matches of the adenine substructure in the input molecule.
    matches = mol.GetSubstructMatches(adenine_mol)
    if not matches:
        return False, "6-aminopurine substructure not found"
    
    # For each match, extract the fragment and compare it (in canonical form)
    # to the reference adenine_smiles.
    for match in matches:
        # Extract the fragment corresponding to the matched atoms.
        # MolFragmentToSmiles returns a SMILES string for only those atoms.
        frag_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=match, canonical=True)
        if frag_smiles == adenine_smiles:
            return True, "Contains 6-aminopurine (adenine) substructure"
    
    # If none of the matches produces exactly the adenine fragment then reject.
    return False, "6-aminopurine-like substructure found, but its core differs from unmodified adenine"

# Example usage (for testing):
if __name__ == "__main__":
    # A known 6-aminopurine (adenine) containing compound (Triacanthine)
    test_smiles = "CC(C)=CCn1cnc(N)c2ncnc12"
    result, reason = is_6_aminopurines(test_smiles)
    print(f"Result: {result}\nReason: {reason}")