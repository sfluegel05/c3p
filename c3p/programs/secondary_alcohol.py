"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: Secondary Alcohol
Definition:
    A secondary alcohol is a compound in which a hydroxy group (-OH) is attached to a 
    saturated (sp³) carbon atom that is bonded to exactly two other carbon atoms and one hydrogen.
    
In many molecules (especially natural products or nucleosides) an –OH group that matches this 
pattern may be present but other –OH groups (e.g. primary alcohols) may also be present.
Here we require that if any –OH group is present then ALL alcohol groups in the molecule 
must be “secondary” for the molecule to be classified as a secondary alcohol.
    
This implementation uses RDKit to add explicit hydrogens and then iterates over every oxygen atom.
For each oxygen, if it is only bonded to one atom (as expected for an –OH), and that one neighbor 
is a carbon, then we count this as an alcohol group. Then we examine the carbon:
    - It must be sp³.
    - Among its bonds (besides the –OH oxygen) there must be exactly two carbons and one hydrogen.
If at least one –OH is found and every –OH group in the molecule meets the secondary alcohol 
criteria, we call the molecule a secondary alcohol; otherwise, it isn’t.
    
Note that with this “whole‐molecule” criterion many complex molecules (even if they contain a 
secondary –OH) will not be classified as secondary alcohols.
"""

from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol overall based on its SMILES string.
    Here “secondary alcohol” means that the only alcohol groups present are secondary.
    
    Steps:
      1. Parse the SMILES and add explicit hydrogens.
      2. Identify all –OH groups (oxygen atoms with exactly one neighbor, that neighbor being carbon).
      3. For each –OH, check that the carbon is sp³ and has exactly 3 neighbors
         (besides the –OH oxygen these must be exactly 2 carbons and 1 hydrogen).
      4. If at least one alcohol group is found and every –OH group qualifies as secondary,
         classify the molecule as a secondary alcohol.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a secondary alcohol, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are available.
    mol = Chem.AddHs(mol)
    
    total_alcohol = 0
    secondary_alcohol = 0
    
    # Loop over all oxygen atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # Expect an alcohol oxygen to have exactly one neighbor.
        if atom.GetDegree() != 1:
            continue
        
        # In an –OH, the sole neighbor should be a carbon.
        neighbor = atom.GetNeighbors()[0]
        if neighbor.GetAtomicNum() != 6:
            continue
        
        # Count this as an alcohol group.
        total_alcohol += 1
        
        # Check that the carbon is saturated (sp³).
        if neighbor.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue  # not secondary if the carbon is not sp³
        
        # Examine the other neighbors of the carbon.
        carbon_neighbor_count = 0
        hydrogen_neighbor_count = 0
        other_neighbor_count = 0
        
        for nb in neighbor.GetNeighbors():
            # Skip the –OH oxygen atom.
            if nb.GetIdx() == atom.GetIdx():
                continue
            if nb.GetAtomicNum() == 6:
                carbon_neighbor_count += 1
            elif nb.GetAtomicNum() == 1:
                hydrogen_neighbor_count += 1
            else:
                other_neighbor_count += 1
        
        # For a secondary alcohol the carbon should have exactly 3 other neighbors:
        # two carbons and one hydrogen.
        if carbon_neighbor_count == 2 and hydrogen_neighbor_count == 1 and other_neighbor_count == 0:
            secondary_alcohol += 1

    # Decide the classification:
    if total_alcohol == 0:
        return False, "No identifiable -OH groups attached to sp³ carbons found"
    if total_alcohol == secondary_alcohol:
        return True, f"All {total_alcohol} alcohol group(s) are secondary (each -OH is on a carbon with 2 C and 1 H)"
    else:
        return False, (f"Not all identified alcohol groups are secondary "
                       f"(found {secondary_alcohol} secondary out of {total_alcohol} total alcohol groups)")

# Example test (uncomment to try out)
# test_smiles = "CCCC(O)CC"  # butan-2-ol (secondary alcohol)
# print(is_secondary_alcohol(test_smiles))