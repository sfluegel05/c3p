"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: Oximes of aldehydes (aldoxime).
Definition: Aldoximes are oximes derived from aldehydes having the structure R–CH=N–OH,
where the carbon in the C=N bond (originally from the aldehyde) bears exactly one hydrogen
and exactly two heavy neighbors (one being the nitrogen of the oxime).
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
        (bool, str): Tuple containing True and a success message if an aldoxime group is found
                     with the correct substitution pattern, otherwise False and a reason.
    """
    # Parse the SMILES into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Note: We do not add explicit hydrogens here, so that the SMARTS pattern sees the natural
    # implicit hydrogen counts (e.g. carbon will show H1 if it has exactly one hydrogen).
    
    # Define a SMARTS pattern that requires:
    # - a carbon atom (atomic number 6) with exactly one hydrogen (i.e. CH1)
    # - double-bonded to a nitrogen
    # - nitrogen single-bonded to an oxygen that is part of an -OH group.
    aldoxime_pattern = Chem.MolFromSmarts("[#6;H1]=[N]-[OH]")
    if aldoxime_pattern is None:
        return False, "Error creating SMARTS pattern for aldoxime"
    
    # Look for substructure matches.
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "Aldoxime functional group (R-CH=N-OH) not found"
    
    # For each found match, verify the substitution pattern.
    for match in matches:
        # match order: [0]: carbon, [1]: nitrogen, [2]: oxygen (in -OH)
        c_idx, n_idx, o_idx = match
        c_atom = mol.GetAtomWithIdx(c_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)
        o_atom = mol.GetAtomWithIdx(o_idx)
        
        # Double-check that the carbon has exactly one hydrogen.
        # (GetTotalNumHs returns the sum of implicit and explicit hydrogens.)
        if c_atom.GetTotalNumHs() != 1:
            continue

        # Verify that the carbon (the aldehyde-derived carbon) has exactly two heavy neighbors.
        # Heavy neighbors are those that are not hydrogens.
        heavy_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # In an aldehyde, the carbon should be bonded to two heavy atoms:
        # one is the oxime nitrogen and one is the R group that was originally attached.
        if len(heavy_neighbors) != 2 or n_idx not in [nbr.GetIdx() for nbr in heavy_neighbors]:
            continue
        
        # Check that the oxygen in the -OH group has exactly one hydrogen.
        if o_atom.GetTotalNumHs() != 1:
            continue
        
        # If all criteria are met, we have identified an aldoxime group.
        return True, "Molecule contains an aldoxime functional group (R-CH=N-OH) with correct aldehyde substitution"
    
    # If none of the matches pass the criteria, return false with a specific reason.
    return False, ("Oxime-like group found, but none of the matches have the expected aldehyde substitution "
                   "(carbon must have exactly one hydrogen and two heavy neighbors)")

# End of program.