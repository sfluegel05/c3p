"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: Diglyceride – a glyceride in which any two of the three hydroxy groups of glycerol are acylated.
A proper diglyceride should have two ester (acyl) groups and one free hydroxyl (–OH) group on the glycerol backbone.
Note: This implementation uses heuristic substructure counts and may be approximate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride (diacylglycerol) based on its SMILES string.
    A diglyceride has glycerol in which two of the hydroxy groups have been acylated.
    The method uses the heuristic that the molecule should contain exactly two ester bonds 
    (i.e. acyl groups attached via –O–C(=O)–) and exactly one free hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a diglyceride, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns:
    # Ester pattern: an oxygen singly bonded to a carbonyl carbon.
    ester_smarts = "[OX2][CX3](=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    
    # Free hydroxyl group pattern (an oxygen with a hydrogen attached)
    free_oh_smarts = "[OX2H]"
    free_oh_pattern = Chem.MolFromSmarts(free_oh_smarts)

    # Count ester groups (acyl groups) in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    num_esters = len(ester_matches)
    
    # Count free hydroxyl groups in the molecule
    free_oh_matches = mol.GetSubstructMatches(free_oh_pattern)
    num_free_oh = len(free_oh_matches)
    
    # For a diglyceride the heuristic is that there must be exactly 2 ester groups
    # and one free hydroxyl on the glycerol backbone.
    if num_esters != 2:
        return False, f"Found {num_esters} ester group(s) (expected exactly 2 for a diglyceride)"
    
    if num_free_oh != 1:
        return False, f"Found {num_free_oh} free hydroxyl group(s) (expected exactly 1 in a glycerol backbone of a diglyceride)"

    # Optionally, one may check that the molecular weight is in a reasonable range
    # for a diglyceride. (Typically, diglycerides will have a molecular weight above ~200 Da.)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) appears too low for a diglyceride"
    
    # If all our heuristic tests pass, we classify it as a diglyceride.
    return True, "Contains a glycerol backbone with 2 acyl (ester) groups and 1 free hydroxyl group (diglyceride)"

# Example usage:
if __name__ == "__main__":
    test_smiles = "C([C@@](COC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])O"
    result, reason = is_diglyceride(test_smiles)
    print(result, reason)