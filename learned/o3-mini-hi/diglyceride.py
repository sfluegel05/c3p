"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: Diglyceride – a glyceride in which any two of the three hydroxy groups are acylated.
A genuine (diacyl)glycerol should have a glycerol backbone in which exactly two –OH groups are 
acylated (via an ester bond O–C(=O)[*]) while the remaining hydroxyl group stays free.

This implementation attempts to match two common backbone arrangements:
  • sn-1,2-diacylglycerol: [CH2:1](O[C:2](=O)[*:7]) - [CH:3](O[C:4](=O)[*:8]) - [CH2:5]([O:6])
  • sn-1,3-diacylglycerol: [CH2:1](O[C:2](=O)[*:7]) - [CH:3]([O:6]) - [CH2:5](O[C:4](=O)[*:8])
  
The free hydroxyl oxygen is marked with map number 6 in these SMARTS.
If a valid match is found then the molecule is classified as a diglyceride.
Otherwise, additional diagnostics (counts of ester bonds and free –OH groups) are returned.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride (diacylglycerol) based on its SMILES string.
    
    A diglyceride is defined as glycerol in which exactly two out of three hydroxy groups are acylated.
    This implementation attempts to match two common backbone patterns using mapped SMARTS:
      - sn-1,2-diacylglycerol: [CH2:1](O[C:2](=O)[*:7])-[CH:3](O[C:4](=O)[*:8])-[CH2:5]([O:6])
      - sn-1,3-diacylglycerol: [CH2:1](O[C:2](=O)[*:7])-[CH:3]([O:6])-[CH2:5](O[C:4](=O)[*:8])
    In both patterns, the free hydroxyl oxygen is tagged as map number 6.
    
    Args:
        smiles (str): The input SMILES string.
    
    Returns:
        bool: True if the molecule is classified as a diglyceride, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Reject very low molecular weight compounds as unrealistic for a diglyceride (>200 Da roughly)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a realistic diglyceride."
    
    # Define SMARTS patterns for sn-1,2 and sn-1,3 diacylglycerol backbones.
    pattern_sn12 = Chem.MolFromSmarts("[CH2:1](O[C:2](=O)[*:7])-[CH:3](O[C:4](=O)[*:8])-[CH2:5]([O:6])")
    pattern_sn13 = Chem.MolFromSmarts("[CH2:1](O[C:2](=O)[*:7])-[CH:3]([O:6])-[CH2:5](O[C:4](=O)[*:8])")
    
    if pattern_sn12 is None or pattern_sn13 is None:
        return None, None  # Cannot compile SMARTS patterns
    
    # Helper function: given a match tuple and a pattern, look for the atom corresponding to map number 6.
    # Then check that it is an oxygen (atomic num 8) attached only to one heavy atom.
    def is_valid_free_oh(match, query_pattern):
        free_oh_idx = None
        for i in range(query_pattern.GetNumAtoms()):
            if query_pattern.GetAtomWithIdx(i).GetAtomMapNum() == 6:
                free_oh_idx = match[i]
                break
        if free_oh_idx is None:
            return False
        atom = mol.GetAtomWithIdx(free_oh_idx)
        if atom.GetAtomicNum() != 8:
            return False
        # Check that the free -OH oxygen is attached only to one heavy atom.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            return False
        return True
    
    # Try matching both patterns. We disable chirality matching for flexibility.
    matches = []
    matches.extend(mol.GetSubstructMatches(pattern_sn12, useChirality=False))
    matches.extend(mol.GetSubstructMatches(pattern_sn13, useChirality=False))
    
    # If any match passes the free hydroxyl check, classify as a diglyceride.
    for match in matches:
        # We can use either pattern's mapping because both mark the free hydroxyl as map 6.
        if is_valid_free_oh(match, pattern_sn12):
            return True, "Matched diglyceride backbone: two acylated hydroxy groups and one free hydroxyl."
    
    # If no backbone pattern could be matched, prepare a diagnostic by counting ester and free OH groups.
    ester_pattern = Chem.MolFromSmarts("[O;!H0]-C(=O)")
    free_oh_pattern = Chem.MolFromSmarts("[OX2H]")
    num_esters = len(mol.GetSubstructMatches(ester_pattern)) if ester_pattern else 0
    num_free_oh = len(mol.GetSubstructMatches(free_oh_pattern)) if free_oh_pattern else 0
    diag = (f"No matching diglyceride backbone found. Detected {num_esters} ester group(s) and "
            f"{num_free_oh} free hydroxyl group(s).")
    return False, diag

# Example usage for testing:
if __name__ == "__main__":
    test_smiles = "C([C@@](COC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])O"
    result, reason = is_diglyceride(test_smiles)
    print(result, reason)