"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: D-Glucoside 
Definition: Any glucoside in which the glycoside group is derived from D-glucose.
In practice we require that (i) the molecule contains a D-glucopyranose ring 
(typical 6-membered ring with one ring-oxygen, five carbons, and a CH2OH branch 
at one of the ring carbons) in either alpha or beta configuration (unmodified from glucose),
and (ii) that the sugar is linked to an aglycone (i.e. at least one atom of the ring 
bonds to an atom outside the ring). 
Note: This detection heuristic uses SMARTS patterns that cover many common representations,
but some true examples may be missed while a few false positives may still occur.
"""

from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    
    A D-glucoside is defined here as a molecule containing an unmodified D-glucopyranose
    ring (in alpha or beta form) that is connected via a glycosidic bond to an aglycone.
    Because the sugar motif is depicted with some variability in natural products, we use 
    two SMARTS patterns that aim to capture many commonly encountered depictions. 
    
    Improvement over earlier attempt:
      - Two SMARTS patterns (one for beta, one for alpha) are defined.
      - For each match found, we check that at least one atom in the sugar ring bonds 
        to an atom not in the matched substructure (i.e. attached to an aglycone).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a D-glucoside, otherwise False.
        str: A brief reason for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for a D-glucopyranose ring.
    # These patterns roughly capture a 6-membered ring (5 carbons + 1 oxygen) that has the typical
    # substitution pattern of D-glucose: a CH2OH branch at the carbon adjacent to the ring oxygen,
    # plus hydroxyl groups on the other ring carbons.
    #
    # Note: The two distinct patterns are used to catch the common beta and alpha configurations.
    pattern_beta = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O")
    pattern_alpha = Chem.MolFromSmarts("OC[C@@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O")
    
    # Try to match both patterns in the target molecule.
    matches_beta = mol.GetSubstructMatches(pattern_beta)
    matches_alpha = mol.GetSubstructMatches(pattern_alpha)
    all_matches = list(matches_beta) + list(matches_alpha)
    
    if not all_matches:
        return False, "No D-glucopyranose substructure found"
    
    # For each sugar-ring match, verify that it is attached (via at least one bond) to an external group.
    # This is our simple test for a glycosidic bond (as opposed to an isolated sugar molecule).
    for match in all_matches:
        match_set = set(match)
        for atom_idx in match_set:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Go through bonds (neighbors)
            for nbr in atom.GetNeighbors():
                # We ignore hydrogens (not explicitly stored generally in RDKit anyway)
                if nbr.GetIdx() not in match_set:
                    # Found a bond connecting the sugar ring to an atom outside the ring.
                    return True, "Contains a D-glucopyranose ring attached via a glycosidic bond"
                    
    return False, "Found a D-glucopyranose ring but no bond connecting it to an aglycone"

# Example usage (Remove or comment out before production use)
if __name__ == '__main__':
    # Example: beta-D-glucosyl-N-eicosanoylsphingosine
    test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
    result, reason = is_D_glucoside(test_smiles)
    print(result, reason)