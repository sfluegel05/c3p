"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: beta-D-glucosiduronic acid conjugates
A glucosiduronic acid resulting from the formal condensation of any substance with 
beta-D-glucuronic acid to form a glycosidic bond.
"""

from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid conjugate based on its SMILES string.
    It checks for the presence of a beta-D-glucuronic acid moiety that is linked to another moiety (aglycone)
    via a glycosidic bond. If the moiety is present but fully isolated (i.e. the free acid), it classifies it as False.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a beta-D-glucosiduronic acid conjugate, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for beta-D-glucuronic acid moiety.
    # This pattern is based on a beta-D-glucuronic acid where the anomeric oxygen is in beta configuration.
    glucuronic_smarts = "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)C(=O)O"
    glucuronic_fragment = Chem.MolFromSmarts(glucuronic_smarts)
    if glucuronic_fragment is None:
        return False, "Error in constructing SMARTS for beta-D-glucuronic acid"

    # Find substructure matches for the glucuronic acid moiety.
    matches = mol.GetSubstructMatches(glucuronic_fragment)
    if not matches:
        return False, "No beta-D-glucuronic acid moiety found in the molecule"

    # Check if at least one beta-D-glucuronic acid fragment is conjugated:
    # That is, at least one atom in the fragment has a neighbor that is not part of the fragment.
    for match in matches:
        # Convert match (tuple of atom indices) to a set for quick lookup.
        match_set = set(match)
        # Loop over atoms in the matched fragment.
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                # If the neighbor is not in the matched fragment, then the fragment is linked to an aglycone.
                if neighbor.GetIdx() not in match_set:
                    return True, "Found beta-D-glucuronic acid moiety conjugated via a glycosidic bond"

    # If none of the matched fragments are connected to an external group, we consider it as free acid.
    return False, "Beta-D-glucuronic acid moiety found but appears to be a free acid (no glycosidic bond detected)"

# Example usage:
if __name__ == "__main__":
    # Test using one of the provided SMILES (ethyl glucuronide)
    test_smiles = "CCO[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O"
    result, reason = is_beta_D_glucosiduronic_acid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)