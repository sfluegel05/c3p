"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: Alpha-hydroxy ketone
Definition: A ketone containing a hydroxy group on the alpha-carbon relative to the C=O group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined as a ketone (C=O with two carbon neighbors)
    that has at least one hydroxy (OH) group on an alpha carbon (a carbon adjacent to the carbonyl carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the OH groups can be detected.
    mol = Chem.AddHs(mol)
    
    # Define SMARTS pattern for a ketone: a carbonyl carbon (C=O) that is bonded to two carbons.
    # This avoids matching aldehydes (which have at least one hydrogen attached to the carbonyl).
    ketone_pattern = Chem.MolFromSmarts("[#6][C](=O)[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not ketone_matches:
        return False, "No ketone (with carbonyl flanked by carbons) found in the molecule"
    
    # For each ketone group found, check for an alpha carbon with an attached hydroxy group.
    for match in ketone_matches:
        # In the SMARTS pattern "[#6][C](=O)[#6]":
        # match[0] is the carbon attached to the carbonyl carbon,
        # match[1] is the carbonyl carbon,
        # match[2] is the other carbon attached to the carbonyl carbon.
        carbonyl_carbon_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_carbon_idx)
        
        # Iterate over neighbors of the carbonyl carbon.
        # We skip the oxygen involved in the C=O double bond.
        for neighbor in carbonyl_atom.GetNeighbors():
            # We are looking for an alpha carbon, so it must be a carbon atom.
            if neighbor.GetAtomicNum() != 6:
                continue
            alpha_carbon = neighbor
            # Now check if this alpha carbon is bonded to a hydroxy group (an oxygen with a hydrogen).
            for sub_neigh in alpha_carbon.GetNeighbors():
                # Look for an oxygen (atomic num 8)
                if sub_neigh.GetAtomicNum() == 8:
                    # Check if this oxygen carries at least one hydrogen.
                    # With explicit hydrogens added, we can check the count of neighboring H atoms.
                    # (Alternatively, one could search for the SMARTS [OX2H].)
                    if any(h.GetAtomicNum() == 1 for h in sub_neigh.GetNeighbors()):
                        # We found an OH group on an alpha carbon next to a ketone group.
                        return True, "Found a ketone group with an alpha hydroxy substituent"
    # If none of the ketone groups had an alpha OH, then return false.
    return False, "Ketone group(s) found, but no alpha hydroxy group attached to an alpha carbon"
    
# Example usage:
if __name__ == "__main__":
    test_smiles = "OCC([C@H](CC)C)=O"  # (3S)-1-hydroxy-3-methylpentan-2-one
    result, reason = is_alpha_hydroxy_ketone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)