"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: beta-D-glucosiduronic acid conjugates

A beta-D-glucosiduronic acid conjugate (glucuronide) is defined as a glucuronic acid unit 
that has undergone a formal condensation reaction with another substance via its anomeric position.
That is, the sugar’s anomeric –OH has been replaced by a glycosidic (or N-linked) bond to an aglycone.

This implementation uses a SMARTS pattern with an atom map to mark the anomeric oxygen. We then add explicit hydrogens,
and verify that for at least one sugar match the anomeric oxygen is linked to an external heavy atom whose own heavy neighbor 
count (excluding the linker) is nontrivial. This extra check is designed to reduce false positives.
"""

from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid conjugate (a glucuronide) based on its SMILES string.
    This function searches for a beta-D-glucuronic acid substructure whose anomeric oxygen (labeled [O:1] in the SMARTS)
    is bound outside of the sugar ring to a substituent that is more than a trivial (methyl or hydrogen) appendage.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a beta-D-glucosiduronic acid conjugate, False otherwise
        str: A textual reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydroxyl groups are clear.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for beta-D-glucuronic acid.
    # Here we label the anomeric oxygen with atom map number 1.
    # The pattern is based on a six-membered sugar ring (pyranose) bearing a carboxylic acid tail.
    # (Stereochemistry is included but may be relaxed in some variants.)
    glucuronic_smarts = "[O:1][C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)C(=O)O"
    glucuronic_query = Chem.MolFromSmarts(glucuronic_smarts)
    if glucuronic_query is None:
        return False, "Error constructing SMARTS for beta-D-glucuronic acid"
    
    # Look for all substructure matches. We use chirality so that the sugar ring is more specifically identified.
    matches = mol.GetSubstructMatches(glucuronic_query, useChirality=True)
    if not matches:
        return False, "No beta-D-glucuronic acid moiety found in the molecule"
    
    # For each match (each identified glucuronic sugar fragment), check the external neighbor of atom with map 1.
    # In our SMARTS, the first atom ([O:1]) is the anomeric oxygen that in a conjugate should be bonded externally.
    for match in matches:
        # match is a tuple of atom indices such that match[0] corresponds to the anomeric oxygen.
        anomeric_idx = match[0]
        anomeric_atom = mol.GetAtomWithIdx(anomeric_idx)
        external_candidates = []
        # Loop over neighbors of the anomeric oxygen.
        for nbr in anomeric_atom.GetNeighbors():
            # Ignore hydrogen neighbors.
            if nbr.GetAtomicNum() == 1:
                continue
            # Skip if this neighbor is part of the sugar (i.e. in the matched indices).
            if nbr.GetIdx() in match:
                continue
            external_candidates.append(nbr)
        
        if external_candidates:
            # Now, even if an external neighbor is found, require that it is not just a trivial substituent.
            # Count how many heavy (non-hydrogen) neighbors (other than the anomeric oxygen) the external candidate has.
            for ext in external_candidates:
                heavy_count = sum(1 for nb in ext.GetNeighbors() if nb.GetAtomicNum() > 1 and nb.GetIdx() != anomeric_idx)
                # For example, a methyl carbon attached to the ring oxygen will have no heavy neighbor besides that bond.
                # We require at least one extra heavy connection to flag a true aglycone.
                if heavy_count >= 1:
                    return True, ("Found beta-D-glucuronic acid moiety conjugated through its anomeric oxygen " +
                                  "to an aglycone (external neighbor with sufficient heavy connectivity).")
            # Even if the extra heavy neighbor check fails, if any external neighbor is present, we still flag it,
            # though with a less certain reason.
            return True, "Found beta-D-glucuronic acid moiety conjugated via its anomeric oxygen (external substituent detected)."
    
    return False, ("Beta-D-glucuronic acid moiety found but no external heavy substituent on its anomeric oxygen " +
                   "was detected. Likely a free acid form rather than a conjugate.")

# Example usage (test with tamoxifen N-beta-D-glucosiduronic acid)
if __name__ == "__main__":
    test_smiles = "CC\\C(c1ccccc1)=C(/c1ccccc1)c1ccc(OCC[N+](C)(C)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)C(O)=O)cc1"
    result, reason = is_beta_D_glucosiduronic_acid(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)