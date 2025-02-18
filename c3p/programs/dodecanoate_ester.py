"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: any fatty acid ester in which the carboxylic acid component is lauric acid (dodecanoate ester)
Definition: the acyl group is exactly CH3(CH2)10C(=O)O.
We first search for the SMARTS pattern "CCCCCCCCCCCC(=O)O" (which represents 12 carbons in a row with the last 
carbon forming the acyl (CO) group esterified to an oxygen) and then, for each match, we verify that the matched alkyl chain 
is terminal (i.e. it is not embedded in a longer chain).
"""

from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester.
    A dodecanoate (laurate) ester is defined as any fatty acid ester in which 
    the acyl (carboxylic acid) group is lauric acid (CH3(CH2)10C(=O)O).
    This function finds candidate substructure matches for the laurate ester group 
    then verifies that the alkyl chain is terminal rather than just being embedded 
    as a fragment of a longer chain.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule contains a dodecanoate ester (laurate ester) moiety, False otherwise.
        str: Explanation for classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define SMARTS for the laurate ester substructure.
    # The SMARTS "CCCCCCCCCCCC(=O)O" represents the group:
    #   CH3(CH2)10C(=O)O 
    # that is, 12 carbons with the 12th forming the carbonyl-ester moiety.
    smarts = "CCCCCCCCCCCC(=O)O"
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        return False, "Failed to create SMARTS pattern for dodecanoate ester."

    # Get all matches of the SMARTS pattern in the molecule.
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Dodecanoate ester substructure not found."

    # We now filter the matches. Some matches might pick up the dodecanoate fragment
    # as part of a longer alkyl chain. The idea is to require that the first atom in the match 
    # (the terminal CH3 group) has no neighboring carbon outside the match.
    for match in matches:
        # match is a tuple of atom indices corresponding to the query atoms in order.
        # In our SMARTS, atom0 should be the terminal CH3.
        atom0 = mol.GetAtomWithIdx(match[0])
        is_terminal = True
        # Check every neighbor of atom0 that is a carbon.
        for neigh in atom0.GetNeighbors():
            # if the neighbor is a carbon and its index is not in our match,
            # then the terminal CH3 is actually part of a longer chain.
            if neigh.GetAtomicNum() == 6 and neigh.GetIdx() not in match:
                is_terminal = False
                break
        if is_terminal:
            # We found a match where the acyl chain is exactly laurate.
            return True, "Molecule contains a dodecanoate ester (laurate ester) group."
    
    return False, "Dodecanoate ester substructure found but acyl chain is extended beyond laurate."

# Example usage:
if __name__ == "__main__":
    # One test with a known laurate ester: 1-lauroyl-sn-glycerol
    test_smiles = "CCCCCCCCCCCC(=O)OC[C@@H](O)CO"
    result, reason = is_dodecanoate_ester(test_smiles)
    print(result, reason)