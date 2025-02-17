"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: Secondary Alcohol
Definition:
    A secondary alcohol is a compound in which a hydroxy group (-OH) is attached to a 
    saturated (sp³) carbon that is bonded to exactly two other carbon atoms and one hydrogen.
    
    In many molecules an –OH group that matches this pattern may be present along with other –OH groups 
    (for instance primary alcohols, tertiary alcohols or even –OH bound to non-sp³ carbons). Here we 
    require that if any –OH group is present then ALL such –OH groups must be secondary for the molecule to be classified as a secondary alcohol.
    
Our updated approach uses RDKit to add explicit hydrogens. Then for every oxygen atom we identify if it is in an –OH group by checking that it has exactly two neighbors – one hydrogen and one carbon. For the attached carbon we then verify that it is sp³ and that (besides the –OH oxygen) it is connected to exactly two carbon atoms and one hydrogen. Only if at least one –OH is found and every –OH meets this secondary pattern do we return success.
"""

from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    Here “secondary alcohol” means that every hydroxyl group (–OH) in the molecule is attached to a 
    saturated (sp³) carbon that is bonded to exactly two other carbons and one hydrogen.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a secondary alcohol, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that we can count hydrogen atoms.
    mol = Chem.AddHs(mol)
    
    total_OH = 0
    secondary_OH = 0
    
    # Iterate over all atoms looking for oxygen atoms that are part of an OH group.
    for atom in mol.GetAtoms():
        # We only care about oxygen.
        if atom.GetAtomicNum() != 8:
            continue
        
        # For an -OH group (after adding hydrogens), the oxygen should have exactly 2 neighbors:
        # one hydrogen and one other atom.
        if atom.GetDegree() != 2:
            continue
        
        neighbors = atom.GetNeighbors()
        # Identify hydrogen and carbon neighbor.
        h_neighbor = None
        c_neighbor = None
        for nb in neighbors:
            if nb.GetAtomicNum() == 1:
                h_neighbor = nb
            elif nb.GetAtomicNum() == 6:
                c_neighbor = nb
        # If we did not find one hydrogen and one carbon, then skip this oxygen.
        if h_neighbor is None or c_neighbor is None:
            continue
        
        # We have identified an -OH group
        total_OH += 1
        
        # Check the attached carbon: it must be sp³.
        if c_neighbor.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue  # not a secondary alcohol
        
        # For a secondary alcohol, the carbon attached to the -OH should be bound to exactly:
        # 2 other carbons and 1 hydrogen (the -OH makes the 4th substituent).
        carbon_neighbors = c_neighbor.GetNeighbors()
        count_carbon = 0
        count_hydrogen = 0
        # We exclude the oxygen atom (the -OH group) from the counts.
        for nb in carbon_neighbors:
            if nb.GetIdx() == atom.GetIdx():
                continue
            if nb.GetAtomicNum() == 6:
                count_carbon += 1
            elif nb.GetAtomicNum() == 1:
                count_hydrogen += 1
            else:
                # If there are any other heteroatoms attached (other than the OH) then it is not our pattern.
                count_hydrogen += 0  # just ignore; we will not count this as matching.
        # Check if exactly two carbons and one hydrogen are attached.
        if count_carbon == 2 and count_hydrogen == 1:
            secondary_OH += 1
    
    if total_OH == 0:
        return False, "No identifiable -OH groups (in proper -OH form) found"
    
    if total_OH == secondary_OH:
        return True, f"All {total_OH} alcohol group(s) are secondary (each -OH is on a sp³ carbon bound to 2 carbons and 1 hydrogen)"
    else:
        return False, (f"Not all identified alcohol groups are secondary "
                       f"(found {secondary_OH} secondary out of {total_OH} total -OH groups)")

# Example usage:
# test_smiles = "CC(O)CC"  # Butan-2-ol (secondary alcohol)
# print(is_secondary_alcohol(test_smiles))