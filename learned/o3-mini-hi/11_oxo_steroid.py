"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: 11-oxo steroid 
Definition: Any oxo steroid that has an oxo substituent at position 11.
Note: This is a heuristic implementation that checks for the presence of a steroid 
nucleus (approximated by the cyclopentanoperhydrophenanthrene core) and a ketone 
group (C=O) in a ring. Without a full steroid numbering scheme, we assume that if the 
ketone is embedded in the core, it is likely the 11-oxo substituent.
"""

from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 11-oxo steroid based on its SMILES string.
    The molecule must contain a steroid nucleus (four fused rings) and a ketone 
    group (C=O) embedded in that ring system (heuristically interpreted as an 
    oxo substituent in the 11-position).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 11-oxo steroid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS pattern for the steroid nucleus.
    # This pattern approximates the cyclopentanoperhydrophenanthrene core (four fused rings).
    steroid_smarts = "C1CC2CCC3C(C1)CCC4C3CCC2C4"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is None:
        return None, None  # In case pattern creation fails
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found."
    
    # Get one match of the steroid core to use its atom indices.
    steroid_match = mol.GetSubstructMatch(steroid_pattern)
    steroid_atoms = set(steroid_match)
    
    # Define a SMARTS pattern to detect a ketone group that is embedded in a ring.
    # [R] means a ring atom, and [CX3](=O) requires a trigonal (sp2) carbon with a double bond to oxygen.
    ketone_smarts = "[R][CX3](=O)[R]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    if ketone_pattern is None:
        return None, None  # If pattern creation fails
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not ketone_matches:
        return False, "No ketone group embedded in a ring was detected (required for 11-oxo)."
    
    # Check whether any of the ketone groups occur on a carbon that is part of the steroid core.
    # (We assume that if a ring ketone carbon is among the steroid nucleus atoms it is the 11-oxo group.)
    found_11oxo = False
    for match in ketone_matches:
        # In the SMARTS [R][CX3](=O)[R], the ketone carbon is the second atom (index 1 of the match).
        ketone_carbon = match[1]
        if ketone_carbon in steroid_atoms:
            found_11oxo = True
            break

    if not found_11oxo:
        return False, "Ketone group not located on the steroid nucleus (likely not 11-oxo)."
    
    # Heuristically, if a steroid nucleus is present and one of its ring carbons
    # carries a ketone function, we classify the molecule as a 11-oxo steroid.
    return True, "Steroid nucleus found with an embedded ketone group, consistent with a 11-oxo steroid."