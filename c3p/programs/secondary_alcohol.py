"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: Secondary Alcohol
Definition:
    A secondary alcohol is a compound in which a hydroxy group (-OH) is attached to a saturated
    carbon atom (sp³) that is bonded to two other carbon atoms and one hydrogen.
    
This implementation uses RDKit to add explicit hydrogens and then iterates over hydroxyl groups (O with one neighbor).
For each -OH group, it identifies the carbon it’s attached to and checks that:
  - The carbon is sp³.
  - The carbon (after adding explicit hydrogens) has exactly 4 neighbors.
  - Excluding the -OH oxygen, exactly two other neighbors are carbons and one is a hydrogen.
If at least one such match is found, the molecule is classified as a secondary alcohol.
"""

from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    
    A secondary alcohol must have an -OH group attached to a saturated (sp³) carbon that is bonded
    to exactly two carbon atoms and one hydrogen (i.e. the tetrahedral carbon has 4 bonds: one -OH, two C, one H).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule contains at least one secondary alcohol group, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts will be available
    mol = Chem.AddHs(mol)
    
    # Iterate over all atoms in the molecule; we want those that represent a hydroxyl group.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue  # skip non-oxygen atoms
        # We are interested in -OH groups so the oxygen should have only one bond.
        if atom.GetDegree() != 1:
            continue
        # Get the neighbor of the oxygen (should be a carbon)
        neighbor = atom.GetNeighbors()[0]
        if neighbor.GetAtomicNum() != 6:
            continue  # not bound to carbon; not an alcohol group
        # Check that the carbon is saturated (sp3). (Note: explicit Hs have been added.)
        if neighbor.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # For a true secondary alcohol the carbon should have exactly 4 bonds
        # (one to the -OH group, two to other carbon groups, and one to a hydrogen).
        if neighbor.GetDegree() != 4:
            continue
        
        # Count the number of carbon and hydrogen neighbors (excluding the -OH oxygen).
        carbon_count = 0
        hydrogen_count = 0
        for nb in neighbor.GetNeighbors():
            # Skip the oxygen (the -OH group) for this count.
            if nb.GetIdx() == atom.GetIdx():
                continue
            if nb.GetAtomicNum() == 6:
                carbon_count += 1
            elif nb.GetAtomicNum() == 1:
                hydrogen_count += 1
            else:
                # if there is any other substituent, then this carbon is not a typical center of a secondary alcohol.
                break
        else:
            # We were able to count without an unexpected substituent.
            if carbon_count == 2 and hydrogen_count == 1:
                return True, "Contains at least one secondary alcohol group"
    return False, "No secondary alcohol functional group found"

# Example test (uncomment to try out)
# test_smiles = "CC(C)(C)NCC(O)COc1cccc2NC(=O)CCc12"
# print(is_secondary_alcohol(test_smiles))