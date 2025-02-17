"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: Oximes of aldehydes (aldoxime).
Definition: Aldoximes are oximes derived from aldehydes having the structure R–CH=N–OH,
where the carbon of the C=N bond (originally forming the aldehyde) bears exactly one hydrogen.
This program identifies an oxime‐like group and then further verifies that the carbon
attached to the =N–OH has exactly one hydrogen and exactly two heavy (non‐hydrogen) neighbors.
"""

from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime (an oxime of an aldehyde) based on its SMILES string.
    
    The function first adds explicit hydrogens and looks for oxime-like groups using a loose
    SMARTS pattern. Then for each candidate match it checks:
      - The carbon (C) atom bonded via a double bond to the nitrogen has exactly one hydrogen.
      - That carbon is connected to exactly two heavy atoms (its neighbors excluding hydrogens),
        one of which is the oxime nitrogen.
      - The oxygen bonded to that nitrogen (–OH) has exactly one hydrogen.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple with True and a correct classification reason if
            the molecule is an aldoxime, otherwise False and a reason for the failure.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are explicit.
    mol = Chem.AddHs(mol)
    
    # Define a loose SMARTS pattern to identify an oxime moiety.
    # The pattern "[C]=[N]-[O]" catches a carbon double-bonded to a nitrogen that is single-bonded to an oxygen.
    aldoxime_pattern = Chem.MolFromSmarts("[C]=[N]-[O]")
    if aldoxime_pattern is None:
        return False, "Error creating SMARTS pattern for aldoxime"
    
    # Find all substructure matches for the pattern.
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "Aldoxime functional group (R-CH=N-OH) not found"
    
    # Iterate over each found match and perform stricter checks.
    for match in matches:
        # match ordering as defined by the SMARTS:
        # match[0]: the carbon, match[1]: the nitrogen, match[2]: the oxygen.
        c_idx, n_idx, o_idx = match
        c_atom = mol.GetAtomWithIdx(c_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)
        o_atom = mol.GetAtomWithIdx(o_idx)
        
        # Check that the carbon (C) atom has exactly one hydrogen.
        if c_atom.GetTotalNumHs() != 1:
            continue  # Not consistent with an aldehyde-derived carbon.
        
        # Check that the carbon atom has exactly two heavy neighbors (non-hydrogen atoms).
        heavy_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # One heavy neighbor must be the oxime nitrogen. In an aldehyde, there should be exactly 2 heavy neighbors.
        if len(heavy_neighbors) != 2 or not any(nbr.GetIdx() == n_idx for nbr in heavy_neighbors):
            continue
        
        # Check that the oxygen has exactly one hydrogen (making it an OH group).
        if o_atom.GetTotalNumHs() != 1:
            continue
        
        # Passed all criteria: We have identified an aldoxime group.
        return True, ("Molecule contains an aldoxime functional group (R-CH=N-OH) "
                      "with the correct aldehyde substitution pattern")
    
    # If no candidate passed the additional substitution criteria.
    return False, ("Oxime-like group found, but none of the matches have the expected aldehyde substitution "
                   "(carbon must have exactly one hydrogen and two heavy neighbors)")

# End of program.