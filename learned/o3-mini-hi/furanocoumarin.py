"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: furanocoumarin
Definition: Any furochromene that consists of a furan ring fused with a coumarin.
The fusion may occur in different ways in several isomers.
The algorithm searches for a coumarin substructure (benzopyrone) and a furan substructure,
and then checks if they share at least two atoms (i.e. are fused).
"""

from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    
    A furanocoumarin is defined as a furochromene, i.e. a coumarin (benzopyrone) fused with a furan ring.
    The algorithm proceeds in two steps:
      1. Identify the coumarin moiety using a SMARTS pattern that captures its benzopyrone core.
      2. Identify aromatic furan rings using a SMARTS pattern.
    Then, for each coumarin match and each furan match, it checks whether they share at least 2 atom indices,
    which is a heuristic for fusion.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a furanocoumarin, False otherwise.
        str: A textual explanation of the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS for coumarin moiety.
    # Coumarin (2H-1-benzopyran-2-one) can be represented as O=c2oc1ccccc1c2.
    coumarin_smarts = "O=c2oc1ccccc1c2"
    coumarin_pattern = Chem.MolFromSmarts(coumarin_smarts)
    if coumarin_pattern is None:
        return False, "Invalid coumarin SMARTS"
    
    # Define a SMARTS for furan ring.
    # Furan is an aromatic 5-membered ring containing one oxygen: c1occc1.
    furan_smarts = "c1occc1"
    furan_pattern = Chem.MolFromSmarts(furan_smarts)
    if furan_pattern is None:
        return False, "Invalid furan SMARTS"
    
    # Find all matches for the coumarin core.
    coumarin_matches = mol.GetSubstructMatches(coumarin_pattern)
    if not coumarin_matches:
        return False, "No coumarin (benzopyrone) moiety found"
    
    # Find all matches for aromatic furan rings.
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "No furan ring found"
    
    # Check if any furan ring is fused to the coumarin core.
    # Fusion is defined here as sharing at least two atoms.
    for coup_match in coumarin_matches:
        coup_set = set(coup_match)
        for furan_match in furan_matches:
            furan_set = set(furan_match)
            common_atoms = coup_set.intersection(furan_set)
            if len(common_atoms) >= 2:
                return True, "Furan ring is fused with the coumarin core, indicating a furanocoumarin structure"
    
    return False, "No fused furan and coumarin moieties found"