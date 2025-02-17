"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: Diglyceride – a glyceride in which any two of the three hydroxy groups of glycerol are acylated.
A proper diglyceride should have a glycerol backbone with exactly two acyl (ester) groups attached and one free hydroxyl.
We attempt to detect two possible glycerol backbone arrangements:
 - sn-1,2-diacylglycerol: CH2(OC(=O)*) - CH(OC(=O)*) - CH2(O)
 - sn-1,3-diacylglycerol: CH2(OC(=O)*) - CH(O) - CH2(OC(=O)*)
The SMARTS patterns below use "*" as a wildcard in place of "[*]" (which previously led to failed substructure queries).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride (diacylglycerol) based on its SMILES string.
    A diglyceride is defined as glycerol in which any two of the three hydroxy groups are acylated.
    
    This implementation checks for two common backbone patterns:
      - sn-1,2-diacylglycerol: CH2(OC(=O)*) - CH(OC(=O)*) - CH2(O)
      - sn-1,3-diacylglycerol: CH2(OC(=O)*) - CH(O) - CH2(OC(=O)*)
    
    If either pattern is found, along with a reasonable molecular weight, the molecule is classified as a diglyceride.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a diglyceride, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Define SMARTS patterns for the diglyceride backbones.
    # Use "*" as a wildcard rather than "[*]" to ensure the SMARTS compiles.
    pattern_sn12 = Chem.MolFromSmarts("[CH2](OC(=O)*)-[CH](OC(=O)*)-[CH2](O)")
    pattern_sn13 = Chem.MolFromSmarts("[CH2](OC(=O)*)-[CH](O)-[CH2](OC(=O)*)")
    
    # Check that the SMARTS patterns compiled successfully.
    if pattern_sn12 is None or pattern_sn13 is None:
        return None, None  # Cannot proceed if pattern compilation failed.
    
    # Search for matches, ignoring chirality for flexibility.
    match_sn12 = mol.GetSubstructMatches(pattern_sn12, useChirality=False)
    match_sn13 = mol.GetSubstructMatches(pattern_sn13, useChirality=False)
    
    if match_sn12 or match_sn13:
        # Optionally, check the molecular weight to support the classification.
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 200:
            return False, f"Molecular weight ({mol_wt:.1f} Da) appears too low for a diglyceride."
        return True, "Contains a glycerol backbone with exactly 2 acyl (ester) groups and 1 free hydroxyl group."
    
    # For further diagnostics, count general ester groups and free hydroxyl groups.
    ester_pattern = Chem.MolFromSmarts("O[C]=O")
    free_oh_pattern = Chem.MolFromSmarts("[OX2H]")
    
    num_esters = len(mol.GetSubstructMatches(ester_pattern)) if ester_pattern is not None else 0
    num_free_oh = len(mol.GetSubstructMatches(free_oh_pattern)) if free_oh_pattern is not None else 0
    
    diagnostic = (f"No matching diglyceride backbone found. Detected {num_esters} ester group(s) and "
                  f"{num_free_oh} free hydroxyl group(s).")
    return False, diagnostic

# Example usage:
if __name__ == "__main__":
    # Test with a known diglyceride example.
    test_smiles = "C([C@@](COC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])O"
    result, reason = is_diglyceride(test_smiles)
    print(result, reason)