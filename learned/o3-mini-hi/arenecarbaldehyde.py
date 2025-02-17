"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
For example, piperonal ([H]C(=O)c1ccc2OCOc2c1), salicylaldehyde ([H]C(=O)c1ccccc1O), etc.
"""

from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is an aldehyde where the carbonyl carbon is bonded to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for the classification
    """
    
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Flags to indicate if we find any aldehyde group (with C=O) or
    # an aldehyde group attached to an aromatic atom.
    found_aldehyde = False
    found_arene_aldehyde = False
    
    # Iterate over all atoms in the molecule.
    for atom in mol.GetAtoms():
        # We are looking for a carbon atom that is part of an aldehyde group.
        if atom.GetAtomicNum() == 6:
            # Get heavy-atom neighbors (exclude hydrogens).
            heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
            # In an aldehyde functional group, the carbonyl carbon (CHO) usually has
            # exactly two heavy neighbors: one oxygen (double-bonded) and one other (R-).
            if len(heavy_neighbors) == 2:
                oxygen_found = False
                aromatic_neighbor_found = False
                # Iterate over the heavy neighbors and check bond types.
                for nbr in heavy_neighbors:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    # Check for double bond to oxygen:
                    if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                        oxygen_found = True
                    # Check if the other neighbor is part of an aromatic ring.
                    elif nbr.GetIsAromatic():
                        aromatic_neighbor_found = True
                        
                if oxygen_found:
                    # We have an aldehyde-like group (CHO) because the carbon has a double bonded oxygen.
                    found_aldehyde = True
                    # Check if the carbonyl carbon is attached to an aromatic atom.
                    if aromatic_neighbor_found:
                        found_arene_aldehyde = True
                        # We can return True immediately once we detect a qualifying arenecarbaldehyde.
                        return True, "Arenecarbaldehyde functional group detected (aldehyde connected to aromatic moiety)."
    
    # Return appropriate message based on what we found.
    if found_aldehyde and not found_arene_aldehyde:
        return False, "Aldehyde group present, but not attached to an aromatic moiety."
    else:
        return False, "No aldehyde group detected."

# Example usage (uncomment to test):
# print(is_arenecarbaldehyde("[H]C(=O)c1ccc2OCOc2c1"))  # piperonal should return True
# print(is_arenecarbaldehyde("CC(=O)C"))  # An aliphatic aldehyde should return False