"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: furanocoumarin
Definition: Any furochromene formed by a furan ring fused with a coumarin (benzopyrone) core.
The fusion may occur in different isomeric forms.
This algorithm uses multiple SMARTS patterns for the coumarin part to capture potentially substituted or stereochemically complex cores.
"""

from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    
    A furanocoumarin is defined as a furochromene:
      a coumarin (benzopyrone) fused with a furan ring.
      
    The algorithm uses several SMARTS patterns for the coumarin moiety:
      Pattern 1: "c1ccc2C(=O)Oc2c1"                (a classical coumarin)
      Pattern 2: "O=C1OC=CC2=CC=CC=C12"               (capturing a lactone motif)
      Pattern 3: "O=C1C=CC2=CC=CC=C12"                (alternative depiction of benzopyrone)
    And a SMARTS pattern for a furan ring:
      "c1occc1"  (a 5-membered aromatic ring with one oxygen and four carbons)
      
    Fusion is defined as the coumarin and furan fragments sharing at least two common atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a furanocoumarin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the input molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define various SMARTS patterns for the coumarin (benzopyrone) moiety.
    coumarin_smarts_list = [
        "c1ccc2C(=O)Oc2c1",      # classical representation of coumarin
        "O=C1OC=CC2=CC=CC=C12",   # lactone variant
        "O=C1C=CC2=CC=CC=C12"     # alternative benzopyrone depiction
    ]
    
    coumarin_matches = []
    for smarts in coumarin_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue   # skip if SMARTS is invalid
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            coumarin_matches.extend(matches)
    if not coumarin_matches:
        return False, "No coumarin (benzopyrone) moiety found"
    
    # Define the SMARTS pattern for a furan ring (5-membered aromatic ring with one oxygen).
    # This pattern ensures one oxygen and exactly four carbons.
    furan_smarts = "c1occc1"
    furan_pattern = Chem.MolFromSmarts(furan_smarts)
    if furan_pattern is None:
        return False, "Invalid furan SMARTS pattern"
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "No furan ring found"
    
    # Filter furan matches to confirm the ring composition (1 oxygen and 4 carbons).
    valid_furan_matches = []
    for match in furan_matches:
        symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in match]
        if symbols.count("O") == 1 and symbols.count("C") == 4:
            valid_furan_matches.append(match)
    if not valid_furan_matches:
        return False, "No valid furan ring (5-membered with 1 oxygen and 4 carbons) found"
    
    # Check for fusion: the coumarin and furan substructures must share at least 2 atoms.
    for coup_match in coumarin_matches:
        coup_set = set(coup_match)
        for furan_match in valid_furan_matches:
            furan_set = set(furan_match)
            if len(coup_set.intersection(furan_set)) >= 2:
                return True, "Fused coumarin and furan moieties detected indicating a furanocoumarin structure"
    
    return False, "No fused furan and coumarin moieties found"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with one of the provided SMILES strings:
    test_smiles = "O1C(C(OC(=O)C(CC)C)C2=C1C=CC3=C2OC(=O)C=C3)"
    result, reason = is_furanocoumarin(test_smiles)
    print(result, reason)