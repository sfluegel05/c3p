"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: Diglyceride – a glyceride in which any two of the three hydroxy groups are acylated.
A genuine (diacyl)glycerol should have a glycerol backbone in which exactly two –OH groups have been 
acylated (via an ester bond O–C(=O)[*]) and the remaining hydroxyl remains free.
This implementation attempts to match two common backbone arrangements:
 • sn-1,2-diacylglycerol: CH2(OC(=O)[*]) – CH(OC(=O)[*]) – CH2(O)
 • sn-1,3-diacylglycerol: CH2(OC(=O)[*]) – CH(O) – CH2(OC(=O)[*])
The explicit SMARTS patterns below use atom-mapping to mark:
    1,2,5: the three glycerol carbon atoms,
    2 and 4: the acyl oxygen atoms that are bound to carbonyls,
    7 and 8: the carbonyl carbon (with a wildcard tail) atoms, and
    6: the free hydroxyl oxygen.
If a match is found (with a minimum molecular weight threshold) then the molecule is classified 
as a diglyceride. Otherwise additional diagnostics (counting ester bonds and free –OH groups)
are reported.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride (diacylglycerol) based on its SMILES string.
    
    A diglyceride is defined as glycerol in which any two of the three hydroxy groups are acylated.
    This implementation looks for two common backbone patterns:
      - sn-1,2-diacylglycerol: CH2(OC(=O)[*]) - CH(OC(=O)[*]) - CH2(O)
      - sn-1,3-diacylglycerol: CH2(OC(=O)[*]) - CH(O) - CH2(OC(=O)[*])
    The SMARTS patterns require that the acyl groups appear as -OC(=O)[*] and that the free –OH 
    is not substituted with a non-hydrogen heavy atom (i.e. not further acylated or phosphorylated).
    
    Args:
        smiles (str): The input SMILES string.
    
    Returns:
        bool: True if the molecule is classified as a diglyceride, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Require a minimum molecular weight (a rough filter to avoid tiny fragments)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) appears too low for a realistic diglyceride."
    
    # Define SMARTS for sn-1,2-diacylglycerol.
    # The mapping:
    #   [CH2:1](O[C:2](=O)[*:7]) - [CH:3](O[C:4](=O)[*:8]) - [CH2:5]([O:6])
    # Here atom6 is the free hydroxyl oxygen. The [*] acts as a wildcard for the acyl tail.
    pattern_sn12 = Chem.MolFromSmarts("[CH2:1](O[C:2](=O)[*:7])-[CH:3](O[C:4](=O)[*:8])-[CH2:5]([O:6])")
    
    # Define SMARTS for sn-1,3-diacylglycerol.
    # The mapping:
    #   [CH2:1](O[C:2](=O)[*:7]) - [CH:3]([O:6]) - [CH2:5](O[C:4](=O)[*:8])
    pattern_sn13 = Chem.MolFromSmarts("[CH2:1](O[C:2](=O)[*:7])-[CH:3]([O:6])-[CH2:5](O[C:4](=O)[*:8])")
    
    if pattern_sn12 is None or pattern_sn13 is None:
        return None, None  # Cannot continue if patterns fail to compile.
    
    # Search for a matching substructure using either pattern. We disable chirality matching for flexibility.
    matches_sn12 = mol.GetSubstructMatches(pattern_sn12, useChirality=False)
    matches_sn13 = mol.GetSubstructMatches(pattern_sn13, useChirality=False)
    
    # If a candidate match is found, perform an extra check on the free hydroxyl substituent.
    # In our patterns, atom mapping ":6" is the free OH. We ensure that the corresponding oxygen is
    # not further substituted by a heavy atom (which would indicate phosphorylation or extra acylation).
    def is_valid_match(match, pattern):
        # In our mapped pattern, index 6 (by our mapping) corresponds to the free hydroxyl oxygen.
        # Get the atom index from the match tuple in the order of our mapping.
        free_oh_idx = match[5]  # our mapping: 1->0, 2->1, 3->2, 4->3, 5->4, 6->5, etc.
        free_oh_atom = mol.GetAtomWithIdx(free_oh_idx)
        # We expect it to be an oxygen (atomic num 8) with degree 1 (only attached to the glycerol carbon).
        # (Hydrogens are implicit.) If it has additional heavy-atom neighbors, reject the match.
        heavy_neighbors = [nbr for nbr in free_oh_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            return False
        return True

    valid_match = False
    reason_match = ""
    for match in matches_sn12:
        if is_valid_match(match, pattern_sn12):
            valid_match = True
            reason_match = "Contains sn-1,2 glycerol backbone with two esterified groups and one free hydroxyl."
            break
    if not valid_match:
        for match in matches_sn13:
            if is_valid_match(match, pattern_sn13):
                valid_match = True
                reason_match = "Contains sn-1,3 glycerol backbone with two esterified groups and one free hydroxyl."
                break

    if valid_match:
        return True, reason_match

    # If neither pattern matched, provide a diagnostic.
    # (Count general ester bonds and free hydroxyl groups in the molecule.)
    ester_pattern = Chem.MolFromSmarts("[O;H0]-C(=O)")
    free_oh_pattern = Chem.MolFromSmarts("[OX2H]")
    num_esters = len(mol.GetSubstructMatches(ester_pattern)) if ester_pattern is not None else 0
    num_free_oh = len(mol.GetSubstructMatches(free_oh_pattern)) if free_oh_pattern is not None else 0
    diag = (f"No matching diglyceride backbone found. Detected {num_esters} general ester group(s) and "
            f"{num_free_oh} free hydroxyl group(s).")
    return False, diag

# Example usage (testing):
if __name__ == "__main__":
    # Examples to test. (Here we test one of the provided SMILES for a diglyceride.)
    test_smiles = "C([C@@](COC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])O"
    result, reason = is_diglyceride(test_smiles)
    print(result, reason)