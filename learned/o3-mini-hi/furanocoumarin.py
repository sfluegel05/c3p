"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: furanocoumarin
Definition: Any furochromene that consists of a furan ring fused with a coumarin.
The fusion may occur in different ways in several isomers.
This updated algorithm uses a broader coumarin SMARTS pattern and relaxes the fusion criteria.
"""

from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    
    A furanocoumarin is defined as a furochromene:
      a coumarin (benzopyrone) fused with a furan ring.
    
    The function uses:
      - A broader coumarin SMARTS pattern representing the benzopyrone core.
        (SMARTS: "c1ccc2C(=O)Oc2c1")
      - A furan ring pattern (SMARTS: "c1occc1") representing the 5-membered aromatic ring.
    
    Then it checks whether the substructures are fused (i.e., the two rings share at least two atoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a furanocoumarin, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a broader coumarin SMARTS pattern.
    # This pattern corresponds to a coumarin core: a benzene fused with a pyrone.
    coumarin_smarts = "c1ccc2C(=O)Oc2c1"
    coumarin_pattern = Chem.MolFromSmarts(coumarin_smarts)
    if coumarin_pattern is None:
        return False, "Invalid coumarin SMARTS pattern"
    
    coumarin_matches = mol.GetSubstructMatches(coumarin_pattern)
    if not coumarin_matches:
        return False, "No coumarin (benzopyrone) moiety found"
    
    # Define a SMARTS for a furan ring:
    # furan is a 5-membered aromatic ring with exactly one oxygen and four carbons.
    furan_smarts = "c1occc1"
    furan_pattern = Chem.MolFromSmarts(furan_smarts)
    if furan_pattern is None:
        return False, "Invalid furan SMARTS pattern"
    
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "No furan ring found"
    
    # Further filter furan matches for rings having exactly 1 oxygen and 4 carbons.
    valid_furan_matches = []
    for match in furan_matches:
        symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in match]
        if symbols.count("O") == 1 and symbols.count("C") == 4:
            valid_furan_matches.append(match)
    if not valid_furan_matches:
        return False, "No valid furan ring (5-membered with 1 oxygen and 4 carbons) found"
    
    # Check for fusion between the coumarin and furan moieties using a relaxed criterion:
    # a fused system should share at least two atoms.
    for coup_match in coumarin_matches:
        coup_set = set(coup_match)
        for fur_match in valid_furan_matches:
            fur_set = set(fur_match)
            common_atoms = coup_set.intersection(fur_set)
            if len(common_atoms) >= 2:
                return True, "Fused coumarin and furan moieties detected indicating a furanocoumarin structure"
    
    return False, "No fused furan and coumarin moieties found"

# Example usage:
if __name__ == "__main__":
    # Test with one example; additional tests can be added.
    test_smiles = "O1C(C(OC(=O)C(CC)C)C2=C1C=CC3=C2OC(=O)C=C3)"  # Edulisin III example
    result, reason = is_furanocoumarin(test_smiles)
    print(result, reason)