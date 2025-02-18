"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: D-Glucoside (Any glucoside in which the glycoside group is derived from D-glucose)
We look for a D-glucopyranose substructure (as found in many natural glycosides)
and then verify that it is attached to an aglycone (not the whole molecule).
"""

from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is defined as any glucoside in which the glycoside group is derived
    from D-glucose. In practice, we require the molecule to contain a D-glucopyranose
    ring (in either alpha or beta configuration) that is connected to additional moiety.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a D-glucoside, False otherwise.
        str: Reason for classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for a D-glucopyranose ring.
    # These patterns match common depictions of the sugar ring in its beta and alpha forms.
    # Note: The sugar ring is a 6-membered ring with 5 carbons and 1 oxygen,
    # with a CH2OH branch on the carbon next to the ring oxygen.
    # The following SMARTS are examples and may not be perfect but work for many cases.
    pattern_beta = Chem.MolFromSmarts("O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    pattern_alpha = Chem.MolFromSmarts("O[C@H]1O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    
    # Get matches for both patterns
    matches_beta = mol.GetSubstructMatches(pattern_beta)
    matches_alpha = mol.GetSubstructMatches(pattern_alpha)
    all_matches = list(matches_beta) + list(matches_alpha)
    
    if not all_matches:
        return False, "No D-glucopyranose substructure found"
    
    # For each sugar-ring match, verify that it is connected to the rest of the molecule.
    # A proper glycoside will have at least one bond from an atom inside the sugar to an atom outside.
    for match in all_matches:
        match_set = set(match)
        # For every atom in the match, check bonds to atoms not in the match.
        for atom_idx in match_set:
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in match_set:
                    # Found a connecting bond from the sugar to an aglycone.
                    return True, "Contains a D-glucopyranose ring attached via a glycosidic bond"
    
    # If no atom in any matched sugar has a bond to an outside atom then the molecule may be free glucose.
    return False, "Found an isolated D-glucopyranose ring but no glycosidic bond to an aglycone"
    
# Example usage (these lines can be removed in production):
if __name__ == '__main__':
    # Example: beta-D-glucosyl-N-eicosanoylsphingosine
    test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
    result, reason = is_D_glucoside(test_smiles)
    print(result, reason)