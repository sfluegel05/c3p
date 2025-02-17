"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
#!/usr/bin/env python
"""
Classifies: Alpha-hydroxy ketone
Definition: A ketone containing a hydroxy group on the alpha‐carbon relative to the C=O group.
Improvements:
 - Adds explicit hydrogens so the –OH groups are visible.
 - Computes a global oxygen-to-carbon ratio. If unusually high (e.g. >1.0), we assume the molecule
   is not a typical alpha‐hydroxy ketone (often sugar-like or highly oxygenated species).
 - Instead of one monolithic SMARTS query, iterates over candidate ketone carbons (a carbon double bonded
   to oxygen) requiring that this carbon bears exactly two carbon neighbors (i.e. is a ketone and not an aldehyde),
   then checks that at least one adjacent (alpha) carbon carries a hydroxy (–OH) group.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined as a ketone (a carbonyl, C=O, with two carbon substituents)
    that has at least one hydroxy (–OH) group attached to one of the alpha-carbons (the carbons adjacent
    to the carbonyl center).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an alpha-hydroxy ketone, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so –OH groups are visible.
    mol = Chem.AddHs(mol)
    
    # Global filter: count carbons and oxygens.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    # If the O/C ratio is high, assume the molecule is likely a sugar-like or highly polyoxygenated compound.
    if c_count > 0 and (o_count / c_count) > 1.0:
        return False, ("Molecule has a high oxygen/carbon ratio ({} O for {} C); "
                       "likely a sugar-like or polyoxygenated species, not a typical alpha‐hydroxy ketone."
                       .format(o_count, c_count))
    
    # Iterate over all atoms looking for candidate ketone carbons.
    for atom in mol.GetAtoms():
        # Process only carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue
        
        # Check if the carbon is double bonded to oxygen (i.e. is a carbonyl center).
        carbonyl_oxygen = None
        for bond in atom.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    carbonyl_oxygen = nbr
                    break
        if carbonyl_oxygen is None:
            continue  # Not a ketone carbon since no C=O found.
        
        # Check the other substituents on this carbon; for a ketone (and not an aldehyde)
        # the carbonyl carbon must have exactly two carbon substituents (ignoring the oxygen).
        alpha_carbons = []
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            # Skip the carbonyl oxygen.
            if nbr.GetIdx() == carbonyl_oxygen.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                alpha_carbons.append(nbr)
        if len(alpha_carbons) != 2:
            continue  # Not a typical ketone.
        
        # Now for each alpha carbon candidate, check for an -OH substituent.
        for alpha in alpha_carbons:
            for bond in alpha.GetBonds():
                # Consider only single bonds.
                if bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                nbr = bond.GetOtherAtom(alpha)
                # Look for an oxygen.
                if nbr.GetAtomicNum() == 8:
                    # Count explicit hydrogens attached to the oxygen.
                    num_H = sum(1 for nb in nbr.GetNeighbors() if nb.GetAtomicNum() == 1)
                    # If at least one hydrogen is attached, we likely have an -OH.
                    if num_H >= 1:
                        return True, ("Found ketone group (C=O with two carbon substituents) with an alpha "
                                      "hydroxy substituent on one adjacent carbon.")
    # If no candidate substructure found, return False.
    return False, "No alpha hydroxy ketone substructure found"

# Example usage:
if __name__ == "__main__":
    # Test with one known positive: (3S)-1-hydroxy-3-methylpentan-2-one
    test_smiles = "OCC([C@H](CC)C)=O"
    result, reason = is_alpha_hydroxy_ketone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)