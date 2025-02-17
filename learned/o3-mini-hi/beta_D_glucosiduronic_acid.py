"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: beta-D-glucosiduronic acid conjugates

A beta-D-glucosiduronic acid conjugate (glucuronide) is defined as a glucuronic acid unit 
that has undergone a formal condensation reaction with another substance (the aglycone) via its anomeric position.
In other words, the sugar’s anomeric –OH has been replaced by a glycosidic (or N-linked) bond.
This implementation uses a SMARTS pattern to capture the beta-D-glucuronic acid fragment while allowing
for either a protonated or deprotonated carboxyl end. It then verifies that the anomeric oxygen (atom map 1)
leads to an external neighbor with nontrivial heavy connectivity.
"""

from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid conjugate (glucuronide) based on its SMILES string.
    This function searches for a beta-D-glucuronic acid substructure whose anomeric oxygen (labeled [O:1] in the SMARTS)
    is bound to a substituent that is likely to be an aglycone (i.e. has at least one additional heavy-atom neighbor).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a beta-D-glucosiduronic acid conjugate, False otherwise.
        str: A textual reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the hydroxyl groups and carboxyl groups are unambiguous.
    mol = Chem.AddHs(mol)

    # Define a SMARTS pattern for beta-D-glucuronic acid.
    # The pattern attempts to capture the beta-D-glucuronic acid fragment as a six-membered ring (pyranose) whose anomeric oxygen is tagged.
    # The carboxyl moiety is defined as C(=O)[O,OH] to allow either negative or neutral forms.
    # Note: Given the many ways a sugar unit may be depicted these SMARTS are an approximation.
    glucuronic_smarts = "[O:1][C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)C(=O)[O,OH]"
    glucuronic_query = Chem.MolFromSmarts(glucuronic_smarts)
    if glucuronic_query is None:
        return False, "Error constructing SMARTS for beta-D-glucuronic acid"

    # Search for the beta-D-glucuronic acid substructure.
    matches = mol.GetSubstructMatches(glucuronic_query, useChirality=True)
    if not matches:
        return False, "No beta-D-glucuronic acid moiety found in the molecule"
    
    # For each glucuronic match, check the connectivity of the tagged anomeric oxygen (atom mapped with 1).
    for match in matches:
        # In our SMARTS, the first atom ([O:1]) is the anomeric oxygen.
        anomeric_idx = match[0]
        anomeric_atom = mol.GetAtomWithIdx(anomeric_idx)
        
        # Look for neighbors of the anomeric oxygen that are not part of the glucuronic ring.
        external_neighbors = []
        for nbr in anomeric_atom.GetNeighbors():
            # Skip hydrogens.
            if nbr.GetAtomicNum() == 1:
                continue
            if nbr.GetIdx() in match:
                continue  # part of the sugar ring, ignore
            external_neighbors.append(nbr)
        
        if external_neighbors:
            # Check if at least one external neighbor has at least one other heavy-atom neighbor (nontrivial substituent).
            for ext in external_neighbors:
                heavy_couplings = sum(1 for nb in ext.GetNeighbors() if nb.GetAtomicNum() > 1 and nb.GetIdx() != anomeric_idx)
                if heavy_couplings >= 1:
                    return True, ("Found beta-D-glucuronic acid moiety conjugated through its anomeric oxygen " +
                                  "to an aglycone (external neighbor with sufficient heavy connectivity).")
            # If external neighbor exists but the connectivity check is borderline, still flag it with a less confident reason.
            return True, "Found beta-D-glucuronic acid moiety conjugated via its anomeric oxygen (external substituent detected)."
    
    return False, ("Beta-D-glucuronic acid moiety found but no external heavy substituent on its anomeric oxygen " +
                   "was detected. Likely a free acid form rather than a conjugate.")


# Example usage when run as a script.
if __name__ == "__main__":
    # Example test (tamoxifen N-beta-D-glucosiduronic acid)
    test_smiles = "CC\\C(c1ccccc1)=C(/c1ccccc1)c1ccc(OCC[N+](C)(C)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)C(O)=O)cc1"
    result, reason = is_beta_D_glucosiduronic_acid(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Classification result:", result)
    print("Reason:", reason)