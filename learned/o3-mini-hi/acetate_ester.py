"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: Acetate ester — Any carboxylic ester where the acid component is acetic acid.
This module determines if a molecule contains an acetate ester moiety, defined as an ester substructure 
where the acyl part is CH3C(=O)-.
"""

from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is defined as any carboxylic ester where the acid component is acetic acid,
    meaning the molecule must contain an ester group in which the acyl part is CH3C(=O)-.
    
    The function searches for the fragment O-C(=O)-CH3. For each match it checks:
      • The methyl group is truly CH3 (a carbon atom with three hydrogens attached),
      • The oxygen linking the acyl part is not bonded to a phosphorus atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is an acetate ester, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for the acetate ester substructure.
    # We look for an oxygen atom (with degree exactly 2) connected to a carbonyl carbon
    # which is further bonded to a methyl group (CH3).
    acetate_pattern = Chem.MolFromSmarts("[O;D2]-C(=O)-[CH3]")
    
    # Get all substructure matches.
    matches = mol.GetSubstructMatches(acetate_pattern)
    
    if not matches:
        return False, "Does not contain the acetate ester moiety (O-C(=O)-CH3)"
    
    # Loop over each match and check detailed criteria.
    for match in matches:
        # Use only the first three atoms of the match to ensure proper unpacking,
        # even if the returned tuple has extra atoms.
        if len(match) < 3:
            continue  # Skip if match is unexpectedly short
        o_idx, c_idx, m_idx = match[:3]
        
        # Retrieve the oxygen atom (linking oxygen) and the methyl carbon.
        o_atom = mol.GetAtomWithIdx(o_idx)
        m_atom = mol.GetAtomWithIdx(m_idx)
        
        # Check that the methyl atom is a carbon with exactly three attached hydrogens.
        if m_atom.GetAtomicNum() != 6 or m_atom.GetTotalNumHs() != 3:
            continue  # Skip this match if the CH3 condition is not met.
        
        # Optional extra check: ensure that the oxygen atom is not bonded to any phosphorus atoms.
        # This helps avoid phosphoester conflicts.
        interfering = False
        for nb in o_atom.GetNeighbors():
            # Skip the bond to the carbonyl carbon.
            if nb.GetIdx() == c_idx:
                continue
            if nb.GetAtomicNum() == 15:  # Phosphorus atomic number.
                interfering = True
                break
        if interfering:
            continue
        
        # If all checks pass, the molecule qualifies as an acetate ester.
        return True, "Contains the acetate ester moiety (O-C(=O)-CH3)"
    
    # If none of the matches satisfy the detailed criteria, return false.
    return False, "Does not contain a validated acetate ester moiety (O-C(=O)-CH3)"

# Example usage:
if __name__ == "__main__":
    # Testing with methyl acetate as it should return True.
    test_smiles = "COC(C)=O"
    result, reason = is_acetate_ester(test_smiles)
    print(result, reason)