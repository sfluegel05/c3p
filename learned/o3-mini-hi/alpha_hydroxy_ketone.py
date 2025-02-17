"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: Alpha-hydroxy ketone
Definition: A ketone containing a hydroxy group on the alpha‐carbon relative to the C=O group.
Improvement notes:
 - We add explicit hydrogens so that hydroxyl groups (–OH) are visible.
 - Instead of using one complex SMARTS, we iterate over candidate ketone groups.
 - For each carbonyl carbon (a carbon with a double bond to oxygen) that has two carbon substituents
   (to exclude aldehydes), we inspect each adjacent (alpha) carbon.
 - If any alpha–carbon has a single bond to an oxygen that itself is bonded to at least one hydrogen,
   we classify the molecule as an alpha‐hydroxy ketone.
 - Further filtering (e.g. avoiding overly decorated sugar‐like scaffolds) may be added later.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined as a ketone (a carbonyl carbon [C(=O)] with two carbon
    substituents) that has at least one hydroxy (OH) group directly attached to an alpha-carbon
    (i.e. one of the two carbons connected to the carbonyl).
    
    The function parses the molecule, adds explicit hydrogens, and then searches for a candidate 
    carbonyl carbon. It then checks the two alkyl (alpha) substituents to see if one bears an -OH group.

    Args:
        smiles (str): A SMILES string representing the molecule.

    Returns:
        bool: True if the molecule is classified as an alpha-hydroxy ketone, False otherwise.
        str: A message indicating the reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that -OH groups appear explicitly.
    mol = Chem.AddHs(mol)
    
    # Iterate over atoms searching for a carbonyl carbon that qualifies as a ketone.
    for atom in mol.GetAtoms():
        # Look at carbons only.
        if atom.GetAtomicNum() != 6:
            continue

        # Check if the atom has a double bond to an oxygen.
        carbonyl_found = False
        carbonyl_oxygen = None
        for bond in atom.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    carbonyl_found = True
                    carbonyl_oxygen = nbr
                    break
        if not carbonyl_found:
            continue
        
        # To be a ketone (and not an aldehyde), the carbonyl carbon should have two carbon neighbors
        # besides the double-bonded oxygen.
        carbon_neighbors = []
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if carbonyl_oxygen is not None and nbr.GetIdx() == carbonyl_oxygen.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                carbon_neighbors.append(nbr)
        if len(carbon_neighbors) != 2:
            continue  # skip if not a ketone

        # Look over each alpha-carbon neighor.
        for alpha in carbon_neighbors:
            # Check if the alpha-carbon has an attached -OH group.
            for bond in alpha.GetBonds():
                nbr = bond.GetOtherAtom(alpha)
                # We look for an oxygen connected by a single bond.
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.SINGLE:
                    # Check that this oxygen carries at least one hydrogen.
                    # (Using explicit hydrogens added earlier, count hydrogen neighbors.)
                    num_H = sum(1 for sub in nbr.GetNeighbors() if sub.GetAtomicNum() == 1)
                    if num_H >= 1:
                        return True, ("Found ketone group (C=O with two carbon neighbors) with an alpha hydroxy "
                                      "substituent on one adjacent carbon.")
    
    # If we reach here, no qualifying alpha-hydroxy ketone motif was found.
    return False, "No alpha hydroxy ketone substructure found"

# Example usage:
if __name__ == "__main__":
    # Test with one known positive: (3S)-1-hydroxy-3-methylpentan-2-one
    test_smiles = "OCC([C@H](CC)C)=O"
    result, reason = is_alpha_hydroxy_ketone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)