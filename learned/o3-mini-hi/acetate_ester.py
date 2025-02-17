"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: Acetate ester — Any carboxylic ester where the acid component is acetic acid.
An acetate ester must contain a fragment of the form O–C(=O)–CH3.
"""

from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is any carboxylic ester where the acid component is acetic acid,
    meaning that the acyl part is CH3C(=O)-.

    This function looks for the substructure O–C(=O)–CH3. For each match, it checks:
      • That the matching substructure has the proper connectivity,
      • The terminal methyl group is indeed a CH3 (i.e. a carbon atom with three hydrogens),
      • The linking oxygen is not bonded to phosphorus atoms (to avoid misclassification).
    
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
    
    # Define a relaxed SMARTS for the acetate ester substructure.
    # We no longer require the oxygen to have exactly two connections.
    acetate_pattern = Chem.MolFromSmarts("[O]-C(=O)-[CH3]")
    
    # Get all substructure matches.
    matches = mol.GetSubstructMatches(acetate_pattern)
    
    if not matches:
        return False, "Does not contain the acetate ester moiety (O-C(=O)-CH3)"
    
    # Iterate over each match and apply detailed checks.
    for match in matches:
        if len(match) < 3:
            continue  # If match is unexpectedly short, skip it.
        o_idx, c_idx, m_idx = match[:3]
        o_atom = mol.GetAtomWithIdx(o_idx)
        m_atom = mol.GetAtomWithIdx(m_idx)
        
        # Check that the terminal methyl atom really is CH3.
        # (i.e. atomic number 6 and exactly 3 hydrogens attached)
        if m_atom.GetAtomicNum() != 6 or m_atom.GetTotalNumHs() != 3:
            continue
        
        # Optional check: ensure that the linking oxygen is not bonded to a phosphorus atom.
        interfering = False
        for nb in o_atom.GetNeighbors():
            if nb.GetIdx() == c_idx:
                continue  # Skip the carbonyl carbon.
            if nb.GetAtomicNum() == 15:  # Phosphorus atomic number.
                interfering = True
                break
        if interfering:
            continue
        
        # If all criteria are met, return True.
        return True, "Contains the validated acetate ester moiety (O-C(=O)-CH3)"
    
    return False, "Does not contain a validated acetate ester moiety (O-C(=O)-CH3)"

# Example usage:
if __name__ == "__main__":
    # Testing with methyl acetate (should return True).
    test_smiles = "COC(C)=O"
    result, reason = is_acetate_ester(test_smiles)
    print(result, reason)