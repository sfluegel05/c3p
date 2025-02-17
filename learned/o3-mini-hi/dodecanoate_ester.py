"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: any fatty acid ester in which the carboxylic acid component is lauric acid (dodecanoate ester)

Definition: The acyl portion should be exactly laurate:
    Lauric acid (dodecanoic acid) is CH3(CH2)10COOH, so the acyl group in the ester is CH3(CH2)10C(=O)O.
    Notice that the alkyl chain (before the C(=O)) is 11 carbons long.
    
The strategy:
1. Use a SMARTS pattern for the laurate ester group: "CCCCCCCCCCC(=O)O"
   (i.e. 11 consecutive carbons followed by the carbonyl ester group).
2. For each match, check that the first atom (the terminal CH3) is not connected to any extra carbon 
   not included in the match. This ensures that the chain is terminal and not embbedded in a larger 
   chain fragment.
"""

from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule contains a dodecanoate (laurate) ester group.
    A dodecanoate ester here is defined as an ester in which the acyl (carboxylic acid) part is lauric acid,
    having the structure CH3(CH2)10C(=O)O.
    
    The procedure is as follows:
      - We first parse the molecule from its SMILES representation.
      - We then search for the substructure defined by the SMARTS "CCCCCCCCCCC(=O)O".
        (This pattern represents 11 aliphatic carbons (CH3(CH2)9CH2) directly attached to the
         carbonyl group, so that total the acid has 12 carbons.)
      - For each match we verify that the terminal carbon (i.e. the methyl group at the beginning)
        is not connected to any carbon atom outside the match. This step prevents a case where the
        dodecanoate fragment is only a part of a longer carbon chain.
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule contains a dodecanoate ester (laurate ester) group, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS pattern for the laurate ester subgroup.
    # NOTE: Lauric acid (dodecanoic acid) is CH3(CH2)10COOH; hence the acyl part CH3(CH2)10C(=O)O
    # should be represented by 11 carbon atoms followed by the carbonyl-ester fragment.
    smarts = "CCCCCCCCCCC(=O)O"  # 11 C's (terminal chain) + C(=O)O
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        return False, "Failed to create SMARTS pattern for dodecanoate ester."

    # Retrieve all substructure matches.
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Dodecanoate ester substructure not found."

    # For each match, ensure that the match is terminal.
    # We interpret the first atom in our SMARTS (the first "C") as the terminal CH3.
    for match in matches:
        # match is a tuple of atom indices; for our SMARTS, match[0] corresponds to the terminal CH3.
        atom0 = mol.GetAtomWithIdx(match[0])
        is_terminal = True
        for neighbor in atom0.GetNeighbors():
            # If neighbor is a carbon and is not inside the matched indices then the chain continues.
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in match:
                is_terminal = False
                break
        if is_terminal:
            return True, "Molecule contains a dodecanoate ester (laurate ester) group."
    
    return False, "Dodecanoate ester substructure found but acyl chain is extended beyond laurate."

# Example usage:
if __name__ == "__main__":
    # One test: 1-lauroyl-sn-glycerol (a known laurate ester)
    test_smiles = "CCCCCCCCCCCC(=O)OC[C@@H](O)CO"
    result, reasoning = is_dodecanoate_ester(test_smiles)
    print(result, reasoning)