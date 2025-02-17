"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: Diglyceride – a glyceride in which any two of the three hydroxy groups of glycerol are acylated.
A proper diglyceride should have a glycerol backbone with exactly two acyl (ester) groups 
attached via –O–C(=O)– and one free hydroxyl (–OH) group.
This implementation attempts to find a 3-carbon substructure corresponding to a glycerol unit in either the 
sn-1,2-diacylglycerol or sn-1,3-diacylglycerol form.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride (diacylglycerol) based on its SMILES string.
    Diglycerides have a glycerol backbone in which exactly two of the three hydroxyl groups are acylated.
    
    We search for one of two SMARTS patterns:
     - sn-1,2-diacylglycerol: CH2(OC(=O)R) - CH(OC(=O)R) - CH2(OH)
     - sn-1,3-diacylglycerol: CH2(OC(=O)R) - CH(OH) - CH2(OC(=O)R)
    These patterns are approximate but help constrain the search to the glycerol backbone.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a diglyceride, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the two common diglyceride backbones.
    # Pattern for sn-1,2-diacylglycerol: free hydroxyl on a terminal CH2.
    pattern_sn12 = Chem.MolFromSmarts("[CH2](OC(=O)[*])-[CH](OC(=O)[*])-[CH2](OH)")
    # Pattern for sn-1,3-diacylglycerol: free hydroxyl on the middle carbon.
    pattern_sn13 = Chem.MolFromSmarts("[CH2](OC(=O)[*])-[CH](OH)-[CH2](OC(=O)[*])")
    
    match_sn12 = mol.GetSubstructMatches(pattern_sn12)
    match_sn13 = mol.GetSubstructMatches(pattern_sn13)

    if match_sn12 or match_sn13:
        # Optional: check that the molecular weight is in a reasonable range for a diglyceride (>200 Da).
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 200:
            return False, f"Molecular weight ({mol_wt:.1f} Da) appears too low for a diglyceride"
        return True, "Contains a glycerol backbone with 2 acyl (ester) groups and 1 free hydroxyl group (diglyceride)"
    
    # If neither backbone pattern is found, provide a diagnostic message.
    # We can also inspect the counts of ester bonds and free -OH groups (note: these counts may be misleading overall).
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)")
    free_oh_pattern = Chem.MolFromSmarts("[OX2H]")
    num_esters = len(mol.GetSubstructMatches(ester_pattern))
    num_free_oh = len(mol.GetSubstructMatches(free_oh_pattern))
    
    return False, (f"No matching glycerol backbone found. "
                   f"Found {num_esters} ester group(s) and {num_free_oh} free hydroxyl group(s), "
                   "which does not match the expected pattern for a diglyceride.")

# Example usage:
if __name__ == "__main__":
    # A known diglyceride example: DG(18:3(6Z,9Z,12Z)/18:1(11Z)/0:0)
    test_smiles = "C([C@@](COC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])O"
    result, reason = is_diglyceride(test_smiles)
    print(result, reason)