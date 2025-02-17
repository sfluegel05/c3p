"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: Acetate ester — Any carboxylic ester where the acid (acyl) component is acetic acid.

An acetate ester must contain an ester substructure in which the acyl piece is CH3C(=O)-.
That is, the molecule must have an O–C(=O)–C moiety where the acyl carbon is attached to exactly 3 hydrogens.
"""

from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is any carboxylic ester where the acyl part is CH3C(=O)-.
    This function searches for the substructure O–C(=O)–C where the acyl carbon (C) is a methyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a validated acetate ester moiety, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a more relaxed SMARTS pattern to catch acetate ester candidates.
    # Instead of "[O]-C(=O)-[CH3]", use "[O]-C(=O)-[C]" so that implicit hydrogens are not an issue.
    acetate_candidate = Chem.MolFromSmarts("[O]-C(=O)-[C]")
    matches = mol.GetSubstructMatches(acetate_candidate)
    
    if not matches:
        return False, "Does not contain the acetate ester moiety (O-C(=O)-CH3)"
    
    # Iterate over each match and check that the acyl carbon is indeed a methyl group.
    for match in matches:
        # Expect match[0]: oxygen, match[1]: carbonyl carbon, match[2]: acyl carbon (candidate methyl)
        if len(match) < 3:
            continue  # unexpected match, skip
        
        o_idx, carbonyl_idx, acyl_idx = match[:3]
        acyl_atom = mol.GetAtomWithIdx(acyl_idx)
        o_atom = mol.GetAtomWithIdx(o_idx)
        
        # Check that the acyl atom is a carbon.
        if acyl_atom.GetAtomicNum() != 6:
            continue
        
        # For a methyl group, the carbon must have exactly 3 hydrogens.
        # RDKit calculates implicit and explicit hydrogens when using GetTotalNumHs().
        if acyl_atom.GetTotalNumHs() != 3:
            continue
        
        # Optionally, ensure that the linking oxygen is not bonded to phosphorus (or other interfering atoms).
        interfering = False
        for nb in o_atom.GetNeighbors():
            # Skip the carbonyl carbon which is part of the acetate candidate.
            if nb.GetIdx() == carbonyl_idx:
                continue
            # For safety, if the oxygen is bound to phosphorus, skip
            if nb.GetAtomicNum() == 15:
                interfering = True
                break
        if interfering:
            continue
        
        # If the candidate passes all checks, we have a valid acetate ester.
        return True, "Contains the validated acetate ester moiety (O-C(=O)-CH3)"
    
    return False, "Does not contain a validated acetate ester moiety (O-C(=O)-CH3)"


# Example usage:
if __name__ == "__main__":
    # Test with methyl acetate, which should be classified as a valid acetate ester.
    test_smiles = "COC(=O)C"
    result, reason = is_acetate_ester(test_smiles)
    print(result, reason)