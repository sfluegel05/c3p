"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
#!/usr/bin/env python
"""
Classifies: Alpha-hydroxy ketone
Definition: A ketone containing a hydroxy group on the alpha‐carbon relative to the C=O group.
Improvements:
 - Adds explicit hydrogens so the –OH groups are visible.
 - Uses a global oxygen/carbon ratio (allowing a slightly higher threshold of 1.2) 
   to filter out sugar-like or highly polyoxygenated compounds.
 - Instead of one monolithic SMARTS query, iterates over candidate ketone carbons
   (a carbon double bonded to oxygen) requiring that this carbon bears exactly
   two carbon substituents (i.e. is a ketone rather than an aldehyde),
   then checks that at least one adjacent (alpha) carbon carries a hydroxy (-OH)
   group. To be sure it is a true hydroxyl group, the oxygen must be bound by a
   single bond and carry at least one explicit hydrogen.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined as a ketone (a carbonyl, C=O, with two carbon substituents)
    that has at least one hydroxy (-OH) group attached to one of the alpha-carbons (the carbons adjacent
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
    
    # Add explicit hydrogens so that -OH groups appear explicitly.
    mol = Chem.AddHs(mol)
    
    # Global filter: count total carbon and oxygen atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    # Allow up to 1.2 oxygen per carbon; sugar-like compounds usually exceed this.
    if c_count > 0 and (o_count / c_count) > 1.2:
        return False, ("Molecule has a high oxygen/carbon ratio ({} O for {} C); "
                       "likely a sugar-like or highly polyoxygenated species, not a typical alpha‐hydroxy ketone."
                       .format(o_count, c_count))
    
    # Iterate over atoms to find candidate ketone carbons.
    for atom in mol.GetAtoms():
        # Process only carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue
        
        # Look for a double bond to an oxygen to identify a carbonyl.
        carbonyl_oxygen = None
        for bond in atom.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    carbonyl_oxygen = nbr
                    break
        if carbonyl_oxygen is None:
            continue  # No carbonyl found for this carbon.
        
        # Check that this carbonyl carbon has exactly two carbon neighbors (ketone vs aldehyde)
        alpha_atoms = []
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            # Skip the carbonyl oxygen.
            if nbr.GetIdx() == carbonyl_oxygen.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                alpha_atoms.append(nbr)
        if len(alpha_atoms) != 2:
            continue  # Not a typical ketone center.
        
        # For each alpha (adjacent) carbon, look for an –OH substituent.
        for alpha in alpha_atoms:
            # Optionally, we ensure that this alpha carbon is not overly oxygenated
            # (for example, if it has more than one oxygen neighbor aside from the ketone bond).
            o_neighbors = [n for n in alpha.GetNeighbors() if n.GetAtomicNum() == 8]
            if len(o_neighbors) > 1:
                continue

            for bond in alpha.GetBonds():
                # Only consider single bonds for a hydroxy attachment.
                if bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                nbr = bond.GetOtherAtom(alpha)
                if nbr.GetAtomicNum() == 8:
                    # To confirm this is a hydroxy group, check that:
                    # 1) The oxygen is connected via a single bond.
                    # 2) It has at least one explicit hydrogen attached.
                    num_H = sum(1 for nb in nbr.GetNeighbors() if nb.GetAtomicNum() == 1)
                    if num_H >= 1:
                        return True, ("Found ketone group (C=O with two carbon substituents) with an alpha "
                                      "hydroxy substituent on one adjacent carbon.")
    
    return False, "No alpha hydroxy ketone substructure found"

# Example usage:
if __name__ == "__main__":
    # Test with one known positive.
    test_smiles = "OCC([C@H](CC)C)=O"  # (3S)-1-hydroxy-3-methylpentan-2-one
    result, reason = is_alpha_hydroxy_ketone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)