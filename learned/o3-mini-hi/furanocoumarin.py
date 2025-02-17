"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: furanocoumarin
Definition: Any furochromene that consists of a furan ring fused with a coumarin.
The fusion may occur in different ways in several isomers.
This updated algorithm employs multiple SMARTS patterns for the coumarin (benzopyrone) core.
"""

from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    
    A furanocoumarin is defined as a furochromene:
      a coumarin (benzopyrone) fused with a furan ring.
      
    The algorithm uses:
      - Two SMARTS patterns for the coumarin moiety:
         Pattern 1: "O=C1OC=CC2=CC=CC=C12"
         Pattern 2: "c1ccc2C(=O)Oc2c1"
      - A furan ring SMARTS pattern: "c1occc1" (ensuring a 5-membered aromatic ring with one oxygen)
      
    It then checks for fusion between a coumarin and a furan (i.e. that they share at least two atoms)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a furanocoumarin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for coumarin moiety:
    coumarin_smarts_list = [
        "O=C1OC=CC2=CC=CC=C12",  # pattern based on the lactone structure
        "c1ccc2C(=O)Oc2c1"       # more classical coumarin depiction
    ]
    coumarin_matches = []
    for smarts in coumarin_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # skip invalid SMARTS
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            coumarin_matches.extend(matches)
    if not coumarin_matches:
        return False, "No coumarin (benzopyrone) moiety found"
    
    # Define SMARTS pattern for a furan ring:
    # furan: 5-membered aromatic ring with exactly one oxygen and four carbons.
    furan_smarts = "c1occc1"
    furan_pattern = Chem.MolFromSmarts(furan_smarts)
    if furan_pattern is None:
        return False, "Invalid furan SMARTS pattern"
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "No furan ring found"
    
    # Filter furan matches to ensure exactly 1 oxygen and 4 carbons.
    valid_furan_matches = []
    for match in furan_matches:
        symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in match]
        if symbols.count("O") == 1 and symbols.count("C") == 4:
            valid_furan_matches.append(match)
    if not valid_furan_matches:
        return False, "No valid furan ring (5-membered with 1 oxygen and 4 carbons) found"

    # Check for fusion between coumarin and furan moieties.
    # Fusion means that the coumarin ring system and the furan ring share at least two atoms.
    for coup_match in coumarin_matches:
        coup_set = set(coup_match)
        for furan_match in valid_furan_matches:
            furan_set = set(furan_match)
            if len(coup_set.intersection(furan_set)) >= 2:
                return True, "Fused coumarin and furan moieties detected indicating a furanocoumarin structure"
    
    return False, "No fused furan and coumarin moieties found"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test SMILES: Edulisin III (example given among many)
    test_smiles = "O1C(C(OC(=O)C(CC)C)C2=C1C=CC3=C2OC(=O)C=C3)"
    result, reason = is_furanocoumarin(test_smiles)
    print(result, reason)