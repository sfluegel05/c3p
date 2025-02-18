"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: any fatty acid ester in which the carboxylic acid component is lauric acid,
i.e. a dodecanoate (laurate) ester. This specifically means that the acyl chain must be exactly
CH3(CH2)10C(=O)O.

This version uses RDKit’s SMARTS repetition syntax ([CH2;H2]){10} to require exactly 10 CH2 groups
between a terminal methyl group ([CH3;H3]) and an ester carbonyl group. Deviations (such as branching or
chain extension) will result in different hydrogen counts or repetition numbers, preventing a match.
"""

from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule contains a dodecanoate (laurate) ester group.
    A dodecanoate ester is defined as having an acyl moiety that is exactly:
      CH3(CH2)10C(=O)O

    The method proceeds as follows:
      1. The input SMILES is converted into an RDKit molecule.
      2. A SMARTS pattern is used that employs repetition syntax to demand:
           - A terminal methyl group ([CH3;H3]) (degree=1 automatically since extra bonds lower H count)
           - Followed by exactly 10 methylene groups ([CH2;H2]) via {10} repetition.
           - Followed immediately by a carbonyl connected to an oxygen, i.e. C(=O)O.
      3. If this exact substructure is found anywhere in the molecule, it is classified as a laurate ester.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains an exact laurate ester group, False otherwise.
        str: Reason explaining the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # We now use a SMARTS string with repetition syntax:
    #   [CH3;H3]-([CH2;H2]){10}-C(=O)O
    # This pattern requires:
    #   * A terminal CH3 with exactly 3 hydrogen atoms (thus only one bond)
    #   * Exactly ten consecutive CH2 (each with exactly 2 hydrogens)
    #   * Followed by a carbonyl carbon bonded to an oxygen in an ester linkage.
    smarts = "[CH3;H3]-([CH2;H2]){10}-C(=O)O"
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        return False, "Could not build SMARTS pattern for laurate ester."

    # Attempt to match the query to the molecule.
    matches = mol.GetSubstructMatches(query)
    if matches:
        return True, "Molecule contains a dodecanoate (laurate) ester group (exact acyl chain CH3(CH2)10C(=O)O)."
    else:
        return False, "No laurate ester group with an exact CH3(CH2)10C(=O)O acyl chain found."

# Example usage:
if __name__ == "__main__":
    # Test with a known laurate ester: 1-lauroyl-sn-glycerol.
    test_smiles = "CCCCCCCCCCCC(=O)OC[C@@H](O)CO"
    result, reasoning = is_dodecanoate_ester(test_smiles)
    print(result, reasoning)