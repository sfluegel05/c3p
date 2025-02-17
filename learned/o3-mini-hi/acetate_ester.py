"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: Acetate ester — Any carboxylic ester where the carboxylic acid component is acetic acid.
This program checks whether the molecule contains an ester fragment having an acyl part of CH3C(=O)-.
We improve upon the previous version by verifying that the methyl group is truly CH3
and ensuring that the oxygen linking the acyl part is not bonded to interfering elements (e.g. phosphorus).
"""

from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is defined as any carboxylic ester where the acid component is acetic acid,
    meaning the molecule must contain an ester group in which the acyl part is CH3C(=O)-.
    
    This function searches for the fragment O-C(=O)-CH3. Then for each match it verifies:
      • the methyl carbon truly represents a CH3 group,
      • the oxygen linking to the carbonyl is not bonded to a phosphorus atom (which might indicate a phosphoester).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is an acetate ester, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the basic acetate ester substructure:
    # This pattern looks for an oxygen (with two connections) linked to a carbonyl carbon,
    # which in turn is bonded to a methyl group.
    acetate_pattern = Chem.MolFromSmarts("[O;D2]-C(=O)-[CH3]")
    matches = mol.GetSubstructMatches(acetate_pattern)
    
    if not matches:
        return False, "Does not contain the acetate ester moiety (O-C(=O)-CH3)"
    
    # Now check each match for additional characteristics.
    for match in matches:
        # match is a tuple (idx_oxygen, idx_carbonyl, idx_methyl)
        o_idx, c_idx, m_idx = match
        o_atom = mol.GetAtomWithIdx(o_idx)
        m_atom = mol.GetAtomWithIdx(m_idx)
        
        # Check that m_atom is truly a methyl: atomic number 6 and exactly three hydrogens are attached.
        # (GetTotalNumHs returns both implicit and explicit hydrogens.)
        if m_atom.GetAtomicNum() != 6 or m_atom.GetTotalNumHs() != 3:
            continue  # Skip this match if methyl is not exactly CH3.
        
        # Optional extra check: ensure the linking oxygen is not bonded to phosphorus 
        # (common in phosphoester units that we would like to exclude).
        o_neighbors = o_atom.GetNeighbors()
        interfering = False
        for nb in o_neighbors:
            if nb.GetIdx() == c_idx:
                continue  # skip the carbonyl carbon
            if nb.GetAtomicNum() == 15:  # phosphorus
                interfering = True
                break
        if interfering:
            continue
        
        # If we pass all checks, we accept the compound as containing an acetate ester group.
        return True, "Contains the acetate ester moiety (O-C(=O)-CH3)"
    
    # If no match passes our stricter criteria, then we classify the molecule as lacking the acetate ester.
    return False, "Does not contain a validated acetate ester moiety (O-C(=O)-CH3)"

# For testing, you might call:
if __name__ == "__main__":
    # Example: methyl acetate should be True.
    test_smiles = "COC(C)=O"
    result, reason = is_acetate_ester(test_smiles)
    print(result, reason)