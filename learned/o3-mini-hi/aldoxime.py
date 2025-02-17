"""
Classifies: CHEBI:22307 aldoxime
"""
#!/usr/bin/env python
"""
Classifies: Oximes of aldehydes (aldoxime).
Definition: Aldoximes are oximes derived from aldehydes having the structure R–CH=N–OH,
where the carbon of the C=N bond (originally forming the aldehyde) bears exactly one hydrogen.
This program checks not only that an oxime-like group is present but also that the carbon
attaches to only one non-hydrogen substituent (besides the oxime nitrogen), as expected from an aldehyde.
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime (an oxime of an aldehyde) based on its SMILES string.
    
    In an aldoxime the functional group is R–CH=N–OH. After adding explicit hydrogens,
    the aldehyde-derived carbon should be marked as having exactly one hydrogen and exactly two
    heavy neighbors: one is the oxime nitrogen and the other is the R substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains an aldoxime functional group that satisfies the criteria, False otherwise.
        str: A reason explaining the classification.
    """
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are accurate in our checking.
    mol = Chem.AddHs(mol)
    
    # Define the SMARTS pattern for an oxime group.
    # This looks for a carbon with exactly one hydrogen double bonded to a nitrogen,
    # which is single bonded to an oxygen atom that itself has one hydrogen (i.e. –OH).
    aldoxime_pattern = Chem.MolFromSmarts("[C;H1](=[N][O;H1])")
    if aldoxime_pattern is None:
        return False, "Error creating SMARTS pattern for aldoxime"
    
    # Get all substructure matches.
    matches = mol.GetSubstructMatches(aldoxime_pattern)
    if not matches:
        return False, "Aldoxime functional group (R-CH=N-OH) not found"
    
    # For each match, check the substitution pattern of the carbon.
    # The candidate carbon (first atom in match) should have exactly two non-hydrogen neighbors:
    # one is the oxime nitrogen (second atom in match) and one is the R substituent.
    for match in matches:
        c_idx, n_idx, o_idx = match
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Count heavy (non-hydrogen) neighbors for the carbon.
        heavy_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # We require exactly 2 heavy neighbors.
        if len(heavy_neighbors) != 2:
            continue  # Not the correct substitution pattern for an aldehyde-derived carbon.
        # Check that one of these neighbors is the nitrogen in our oxime group.
        if not any(nbr.GetIdx() == n_idx for nbr in heavy_neighbors):
            continue  # Does not have the required nitrogen neighbor.
        # Passed the stricter criteria: carbon is bound to exactly the oxime nitrogen and one R substituent.
        return True, "Molecule contains an aldoxime functional group (R-CH=N-OH) with the correct aldehyde substitution pattern"
    
    # No match passed the additional aldehyde substitution requirement.
    return False, "Aldoxime-like group found, but the carbon does not have the expected aldehyde substitution (exactly one heavy substituent besides N)"
    
# End of program.