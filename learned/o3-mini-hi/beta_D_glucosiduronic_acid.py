"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: beta-D-glucosiduronic acid conjugates

A glucosiduronic acid resulting from the formal condensation of any substance with 
beta-D-glucuronic acid to form a glycosidic (or sometimes N-linked) bond.
This program attempts to discriminate a glucuronide conjugate (with a beta-D-glucuronic acid
moiety linked via its anomeric position to an aglycone) from the free acid form.
"""

from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid conjugate (a glucuronide)
    based on its SMILES string. In a glucuronide the beta-D-glucuronic acid moiety is linked 
    (by an oxygen or nitrogen) to another moiety (the aglycone) – that is, at the anomeric 
    position the oxygen (or rarely, nitrogen) is not simply bound to the sugar ring and a hydrogen.
    
    The function adds explicit hydrogens so that it can inspect the degree of the anomeric atom.
    Note that this test is heuristic and may misclassify some edge cases.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a beta-D-glucosiduronic acid conjugate, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES and add explicit hydrogens (to make sure –OH vs substituted O are distinguished)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for beta-D-glucuronic acid.
    # We assume that in a beta-D-glucuronic acid moiety the anomeric (linkable) oxygen is the first atom.
    # (Stereochemistry is encoded in the SMARTS.)
    glucuronic_smarts = "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)C(=O)O"
    glucuronic_fragment = Chem.MolFromSmarts(glucuronic_smarts)
    if glucuronic_fragment is None:
        return False, "Error constructing SMARTS for beta-D-glucuronic acid"
    
    matches = mol.GetSubstructMatches(glucuronic_fragment)
    if not matches:
        return False, "No beta-D-glucuronic acid moiety found in the molecule"
    
    # For each match, examine the anomeric oxygen (the first atom of our SMARTS match).
    # In the free acid, that oxygen would be connected only to the ring (one heavy neighbor)
    # and its hydrogen(s). In a conjugate the anomeric oxygen is substituted (i.e. has an extra heavy neighbor).
    for match in matches:
        # match is a tuple of atom indices that match our SMARTS;
        # by construction, the first atom in our SMARTS is the potential anomeric oxygen.
        anomeric_idx = match[0]
        anomeric_atom = mol.GetAtomWithIdx(anomeric_idx)
        # Count heavy atom neighbors NOT in the matched fragment.
        external_neighbors = []
        for nbr in anomeric_atom.GetNeighbors():
            # We are only interested in non-hydrogen atoms
            if nbr.GetAtomicNum() == 1:
                continue
            if nbr.GetIdx() not in match:
                external_neighbors.append(nbr)
        # Debug: In a free acid, we expect 0 external heavy neighbors (the oxygen is only attached to the sugar ring).
        # In a conjugate, the linking atom should have at least one external heavy neighbor.
        if external_neighbors:
            # Optionally, we can check that the attached aglycone is not trivially small (e.g. a single methyl)
            # Here, we count heavy atoms in the fragment attached via this bond.
            # (This is a rough check; one may improve it by more sophisticated fragmentation.)
            for ext in external_neighbors:
                frag = Chem.rdmolops.GetMolFrags(Chem.FragmentOnBonds(mol, [mol.GetBondBetweenAtoms(anomeric_idx, ext).GetIdx()]), asMols=True)
                # We choose the fragment that does NOT contain the glucuronic match:
                ext_frag = None
                for f in frag:
                    f_idxs = f.GetProp("_fragAtomMapNumbers").split() if f.HasProp("_fragAtomMapNumbers") else None
                    # Rather than over-engineer this, we simply check the number of heavy atoms in the neighbor.
                    # If the neighbor itself is more than just a single atom then we consider it a true aglycone.
                    if ext.GetDegree() > 1:
                        return True, "Found beta-D-glucuronic acid moiety conjugated via a glycosidic (or N-linked) bond"
            # Even if the attached fragment is small, we count its presence as a potential conjugation.
            return True, "Found beta-D-glucuronic acid moiety conjugated via a glycosidic (or N-linked) bond"
    
    # If none of the glucuronic fragments were found to be connected to an external heavy atom,
    # we classify it as free acid.
    return False, "Beta-D-glucuronic acid moiety found but appears to be a free acid (no conjugated aglycone detected)"

# Example usage (you can test with any SMILES string):
if __name__ == "__main__":
    # Example: tamoxifen N-beta-D-glucosiduronic acid (previously missed)
    test_smiles = "CC\\C(c1ccccc1)=C(/c1ccccc1)c1ccc(OCC[N+](C)(C)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)C(O)=O)cc1"
    result, reason = is_beta_D_glucosiduronic_acid(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)