"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: Oximes of aldehydes (aldoxime).
Definition: Aldoximes are oximes derived from aldehydes having the structure R–CH=N–OH,
where the carbon in the C=N bond (originally from the aldehyde) bears exactly one hydrogen
and exactly two heavy neighbors (one being the oxime nitrogen).
"""

from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime (oxime of an aldehyde) based on its SMILES string.
    An aldoxime has the structure R–CH=N–OH, where the carbon (CH) double-bonded to a nitrogen
    carries exactly one hydrogen and has exactly two heavy neighbors (one being the oxime nitrogen).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): Tuple containing True and a success message if an aldoxime group is found
                     with the correct substitution pattern; otherwise False and a reason.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are clear.
    mol = Chem.AddHs(mol)
    
    # Create a SMARTS pattern which requires:
    # - A carbon atom with exactly one hydrogen (CH1),
    # - double-bonded to a nitrogen, and
    # - the nitrogen is single-bonded to an oxygen that itself carries exactly one hydrogen (-OH).
    # This pattern is "[C;H1]=[N](-[O;H1])"
    aldoxime_pattern = Chem.MolFromSmarts("[C;H1]=[N](-[O;H1])")
    if aldoxime_pattern is None:
        return False, "Error creating SMARTS pattern for aldoxime"
    
    # Search for the aldoxime substructure.
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "Aldoxime functional group (R-CH=N-OH) not found"
    
    # For each found match, verify that the substitution is correct.
    # The expected match order is: [0]: carbon, [1]: nitrogen, [2]: oxygen.
    for match in matches:
        c_idx, n_idx, o_idx = match
        c_atom = mol.GetAtomWithIdx(c_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)
        o_atom = mol.GetAtomWithIdx(o_idx)
        
        # Ensure none of the atoms carry a formal charge. This helps avoid matching
        # species that are not classical aldoximes (e.g. nitro or aci-nitro compounds).
        if c_atom.GetFormalCharge() != 0 or n_atom.GetFormalCharge() != 0 or o_atom.GetFormalCharge() != 0:
            continue
        
        # Check that the carbon atom (aldehyde-derived carbon) indeed has exactly one hydrogen.
        if c_atom.GetTotalNumHs() != 1:
            continue
        
        # Verify that the carbon has exactly two heavy (non-hydrogen) neighbors.
        heavy_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # One heavy neighbor must be the nitrogen of the oxime.
        if len(heavy_neighbors) != 2 or n_idx not in [nbr.GetIdx() for nbr in heavy_neighbors]:
            continue
        
        # Check the oxygen in the -OH has exactly one hydrogen.
        if o_atom.GetTotalNumHs() != 1:
            continue
        
        # If all criteria are met, the molecule has the correct aldoxime group.
        return True, "Molecule contains an aldoxime functional group (R-CH=N-OH) with correct aldehyde substitution"
    
    # If no match satisfies the strict criteria, report the failure.
    return False, ("Oxime-like group found, but none have the expected aldehyde substitution: "
                   "carbon must have exactly one hydrogen and exactly two heavy neighbors, with neutral atoms.")

# End of program.