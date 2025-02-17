"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: Alpha-hydroxy ketone
Definition: A ketone containing a hydroxy group on the alpha‐carbon relative to the C=O group.
Improvements made:
 - We add explicit hydrogens so that -OH groups are visible.
 - We compute a global O-to-C ratio and, if unusually high (as in many sugar-like or polyoxygenated compounds),
   we then assume the molecule is likely not a typical alpha–hydroxy ketone.
 - Instead of a monolithic SMARTS, we iterate over each carbon that has a double bond to oxygen. For each such
   “ketone carbon”, we require it has two carbon neighbors (to distinguish from aldehydes) and then inspect each
   of those α–carbons for a hydroxyl (–OH) group.
 - These refinements address many of the previously seen false positives while trying to retain true positives.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined as a ketone (a carbonyl [C=O] with two carbon substituents)
    that has at least one hydroxy (OH) group attached to one of the alpha-carbons (the carbons connected
    to the carbonyl center).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an alpha-hydroxy ketone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that -OH groups appear explicitly.
    mol = Chem.AddHs(mol)
    
    # Global filtering: count carbons and oxygens.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    # If the molecule is overly oxygenated (e.g. many sugar-like scaffolds), we caution classification.
    if c_count > 0 and (o_count / c_count) > 1.0:
        return False, ("Molecule has a high oxygen/carbon ratio ({} O for {} C); "
                       "likely a sugar-like or polyoxygenated species, not a typical alpha-hydroxy ketone."
                       .format(o_count, c_count))
    
    # Iterate over atoms looking for candidate ketone carbons.
    for atom in mol.GetAtoms():
        # Process only carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue

        # Look for a double bond to an oxygen (carbonyl).
        has_carbonyl = False
        carbonyl_oxygen = None
        for bond in atom.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    has_carbonyl = True
                    carbonyl_oxygen = nbr
                    break
        if not has_carbonyl:
            continue
        
        # To be a ketone (and not an aldehyde) the carbonyl carbon should have exactly two carbon
        # neighbors not including the double-bonded oxygen.
        alpha_candidates = []
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if carbonyl_oxygen is not None and nbr.GetIdx() == carbonyl_oxygen.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                alpha_candidates.append(nbr)
        if len(alpha_candidates) != 2:
            continue   # skip if not a typical ketone
        
        # For each alpha-carbon candidate check if it bears an -OH substituent.
        for alpha in alpha_candidates:
            for bond in alpha.GetBonds():
                # Look only for single bonds.
                if bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                nbr = bond.GetOtherAtom(alpha)
                # Check for an oxygen that would be part of an -OH group.
                if nbr.GetAtomicNum() == 8:
                    # Count explicit hydrogens attached to the oxygen.
                    num_H = sum(1 for sub in nbr.GetNeighbors() if sub.GetAtomicNum() == 1)
                    # Also verify that the bond is not part of a carbonyl (should be single, not double).
                    if num_H >= 1:
                        return True, ("Found ketone group (C=O with two carbon neighbors) with an alpha hydroxy "
                                      "substituent on one adjacent carbon.")
    # If no candidate motif found, then the molecule is not an alpha-hydroxy ketone.
    return False, "No alpha hydroxy ketone substructure found"

# Example usage:
if __name__ == "__main__":
    # Test with one known positive: (3S)-1-hydroxy-3-methylpentan-2-one
    test_smiles = "OCC([C@H](CC)C)=O"
    result, reason = is_alpha_hydroxy_ketone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)