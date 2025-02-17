"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: Oximes of aldehydes (aldoxime).
Definition: Aldoximes are oximes derived from aldehydes having the structure R–CH=N–OH,
where the carbon in the C=N bond (originally from the aldehyde) bears exactly one hydrogen.
This program identifies an aldoxime group using a SMARTS pattern that explicitly enforces the hydrogen count
and then further verifies that the aldehyde-derived carbon has exactly two heavy (non‐hydrogen) neighbors.
"""

from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime (oxime of an aldehyde) based on its SMILES string.
    An aldoxime has the structure R–CH=N–OH, where the carbon (CH) double-bonded to nitrogen has
    exactly one hydrogen and exactly two heavy neighbors (one must be the oxime nitrogen).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple containing True and a success message if an aldoxime group is found,
                     otherwise False and a message for the failure.
    """
    # Parse the SMILES into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are available.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern that explicitly requires:
    # - a carbon atom (atomic number 6) with exactly one hydrogen (CH1)
    # - double-bonded to a nitrogen
    # - nitrogen single-bonded to an oxygen which is part of an -OH group.
    # This pattern will capture R–CH=N–OH where the CH is correctly substituted.
    aldoxime_pattern = Chem.MolFromSmarts("[#6;H1]=[N]-[OH]")
    if aldoxime_pattern is None:
        return False, "Error creating SMARTS pattern for aldoxime"
    
    # Look for substructure matches.
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "Aldoxime functional group (R-CH=N-OH) not found"
    
    # For each found match, further verify the substitution pattern.
    for match in matches:
        # According to our SMARTS, the match order is:
        #   match[0] : the carbon (CH1)
        #   match[1] : the nitrogen
        #   match[2] : the oxygen (of the OH group)
        c_idx, n_idx, o_idx = match
        c_atom = mol.GetAtomWithIdx(c_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)
        o_atom = mol.GetAtomWithIdx(o_idx)
        
        # Double-check that the carbon has exactly one hydrogen (should be true by SMARTS, but we verify).
        if c_atom.GetTotalNumHs() != 1:
            continue
        
        # Verify that the carbon has exactly two heavy neighbors (non-hydrogen atoms).
        heavy_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # In an aldehyde originally, the carbon should have its two heavy neighbors: one is the oxime nitrogen,
        # and the other is belonging to R.
        if len(heavy_neighbors) != 2 or n_idx not in [nbr.GetIdx() for nbr in heavy_neighbors]:
            continue
        
        # Check that the oxygen atom in the -OH group has exactly one hydrogen.
        if o_atom.GetTotalNumHs() != 1:
            continue
        
        # If all criteria are met, we've successfully identified an aldoxime functional group.
        return True, "Molecule contains an aldoxime functional group (R-CH=N-OH) with correct aldehyde substitution"
    
    # If none of the matches passed the criteria, return false with a specific reason.
    return False, ("Oxime-like group found, but none of the matches have the expected aldehyde substitution "
                   "(carbon must have exactly one hydrogen and two heavy neighbors)")

# End of program.