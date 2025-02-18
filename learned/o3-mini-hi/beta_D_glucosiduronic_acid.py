"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: beta-D-glucosiduronic acid conjugates

A glucosiduronic acid resulting from the formal condensation of any substance with 
beta-D-glucuronic acid to form a glycosidic (or N-linked) bond.
This program attempts to discriminate a glucuronide conjugate (with a beta-D-glucuronic acid
moiety linked via its anomeric position to an aglycone) from the free acid form.
"""

from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid conjugate (a glucuronide)
    based on its SMILES string. A glucuronide has a beta-D-glucuronic acid moiety 
    whose anomeric oxygen has an extra heavy neighbor outside the sugar ring.
    
    The function adds explicit hydrogens to help distinguish between hydroxyl groups 
    on the sugar ring and a substituted oxygen.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a beta-D-glucosiduronic acid conjugate, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define the SMARTS for beta-D-glucuronic acid.
    # The SMARTS specifies the sugar ring with its stereochemistry and carboxylic acid.
    glucuronic_smarts = "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)C(=O)O"
    glucuronic_fragment = Chem.MolFromSmarts(glucuronic_smarts)
    if glucuronic_fragment is None:
        return False, "Error constructing SMARTS for beta-D-glucuronic acid"
        
    matches = mol.GetSubstructMatches(glucuronic_fragment)
    if not matches:
        return False, "No beta-D-glucuronic acid moiety found in the molecule"
    
    # Check each glucuronic acid match.
    # In the free acid the anomeric oxygen (first atom in the match) is only bound within the sugar ring.
    # In a conjugate, that oxygen will be linked to an aglycone (i.e. an extra heavy neighbor).
    for match in matches:
        # The first atom in our SMARTS is assumed to be the anomeric oxygen.
        anomeric_idx = match[0]
        anomeric_atom = mol.GetAtomWithIdx(anomeric_idx)
        external_neighbors = []
        for nbr in anomeric_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                continue
            # If the neighbor's index is not part of the matched glucuronic fragment, add it.
            if nbr.GetIdx() not in match:
                external_neighbors.append(nbr)
        if external_neighbors:
            # Optionally check that the external neighbor is more than just a trivial substituent.
            # Here we require that the external neighbor has at least one heavy neighbor (degree > 1).
            for ext in external_neighbors:
                # Retrieve the bond between the anomeric oxygen and the external neighbor.
                bond = mol.GetBondBetweenAtoms(anomeric_idx, ext.GetIdx())
                # If bond exists and the external neighbor is connected to other atoms, consider it as aglycone.
                if ext.GetDegree() > 1:
                    return True, "Found beta-D-glucuronic acid moiety conjugated through its anomeric oxygen to an aglycone."
            # Even if the external neighbor doesn't have a high degree, we count its presence.
            return True, "Found beta-D-glucuronic acid moiety conjugated via its anomeric oxygen."
    
    # If none of the glucuronic fragments have an external heavy neighbor, it is likely a free acid.
    return False, "Beta-D-glucuronic acid moiety found but appears to be in its free acid form (no conjugated aglycone)."

# Example usage (test with a sample SMILES):
if __name__ == "__main__":
    # Example: tamoxifen N-beta-D-glucosiduronic acid
    test_smiles = "CC\\C(c1ccccc1)=C(/c1ccccc1)c1ccc(OCC[N+](C)(C)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)C(O)=O)cc1"
    result, reason = is_beta_D_glucosiduronic_acid(test_smiles)
    print("Test SMILES:", test_smiles)
    print("Result:", result)
    print("Reason:", reason)