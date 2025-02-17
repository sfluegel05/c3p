"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: furanocoumarin
Definition: Any furochromene that consists of a furan ring fused with a coumarin.
The fusion may occur in different ways in several isomers.
The algorithm attempts to:
  1. Identify the coumarin (benzopyrone) core using a slightly broadened SMARTS pattern.
  2. Identify a furan ring by matching a 5‐membered aromatic ring (with exactly 1 oxygen and 4 carbons).
  3. Ensure that the two substructures are fused (i.e. share exactly 2 atoms).
Additional filtering is applied so that the coumarin match only contains carbon and oxygen.
Note: This is a heuristic approach and may mis‐classify borderline cases.
"""

from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.

    A furanocoumarin is defined as a furochromene:
      a coumarin (benzopyrone) fused with a furan ring.
    This function searches for:
      - A coumarin substructure using a relaxed SMARTS (c1ccc2c(c1)oc(=O)c2)
      - A furan ring (c1occc1) that is a 5-membered aromatic ring composed of exactly one oxygen and four carbons.
    If any coumarin match and any furan match share exactly two atoms (i.e. typical for ring-fusion),
    the molecule is considered a furanocoumarin.
    
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
    
    # Define a more general SMARTS for the coumarin (benzopyrone) moiety.
    # This pattern represents a common coumarin skeleton.
    coumarin_smarts = "c1ccc2c(c1)oc(=O)c2"
    coumarin_pattern = Chem.MolFromSmarts(coumarin_smarts)
    if coumarin_pattern is None:
        return False, "Invalid coumarin SMARTS pattern"
    
    # Find all substructure matches for the coumarin core.
    coumarin_matches = mol.GetSubstructMatches(coumarin_pattern)
    if not coumarin_matches:
        return False, "No coumarin (benzopyrone) moiety found"
    
    # Filter coumarin matches to ensure only C and O atoms are present (to avoid heterocycles with N, etc.)
    valid_coumarin_matches = []
    for match in coumarin_matches:
        sub_atoms = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in match]
        # coumarin core should contain only C and O.
        if all(sym in ("C", "O") for sym in sub_atoms):
            valid_coumarin_matches.append(match)
    if not valid_coumarin_matches:
        return False, "Coumarin moiety found but contains atoms other than C and O"
    
    # Define a SMARTS for a furan ring:
    # furan is a 5-membered aromatic ring with exactly one oxygen atom.
    furan_smarts = "c1occc1"
    furan_pattern = Chem.MolFromSmarts(furan_smarts)
    if furan_pattern is None:
        return False, "Invalid furan SMARTS pattern"
    
    # Find all furan substructure matches.
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "No furan ring found"
    
    # Further filter furan_matches to ensure the ring consists of 1 oxygen and 4 carbon atoms.
    valid_furan_matches = []
    for match in furan_matches:
        symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in match]
        if symbols.count("O") == 1 and symbols.count("C") == 4:
            valid_furan_matches.append(match)
    if not valid_furan_matches:
        return False, "No valid furan ring (5-membered ring with 1 oxygen and 4 carbons) found"
    
    # Check for fusion: a furan ring is considered fused to the coumarin if the two matches share exactly two atoms.
    for coup_match in valid_coumarin_matches:
        coup_set = set(coup_match)
        for fur_match in valid_furan_matches:
            fur_set = set(fur_match)
            common_atoms = coup_set.intersection(fur_set)
            # In a typical fused ring system, the two rings share exactly two atoms.
            if len(common_atoms) == 2:
                return True, "Furan ring is fused with the coumarin core, indicating a furanocoumarin structure"
    
    return False, "No fused furan and coumarin moieties found"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples; feel free to add others.
    test_smiles = "O1C(C(OC(=O)C(CC)C)C2=C1C=CC3=C2OC(=O)C=C3)"  # Edulisin III
    result, reason = is_furanocoumarin(test_smiles)
    print(result, reason)