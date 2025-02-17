"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CDP-diacylglycerol
Definition: A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) at the 1- and 2-positions.
Improved heuristics:
  - Checks for a cytidine moiety.
  - Checks for at least two acyl ester bonds.
  - Checks for diphosphate connectivity (two phosphorus atoms bridged by an oxygen).
  - Additionally verifies the presence of a diacylglycerol fragment.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    
    A CDP-diacylglycerol is defined as a CDP-glycerol (i.e. a cytidine diphosphate
    linked to glycerol) with acyl groups (usually fatty acyl chains) at the 1- and 2-positions.
    
    The function applies the following heuristic steps:
      1. Verify the presence of a cytidine moiety.
      2. Check for at least two acyl ester bonds.
      3. Check for diphosphate connectivity (two P atoms bridged by an O).
      4. Verify the presence of a diacylglycerol fragment, i.e. a three-carbon chain 
         in which two positions are esterified with acyl groups and the third is attached (via oxygen)
         to a phosphate group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a CDP-diacylglycerol, False otherwise.
        str: A message describing the reasoning behind the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for a cytidine moiety.
    # This SMARTS pattern captures the cytosine ring found in cytidine.
    cytidine_smarts = "n1ccc(N)nc1=O"
    cytidine_pattern = Chem.MolFromSmarts(cytidine_smarts)
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "Cytidine moiety not found"
    
    # 2. Check for acyl ester bonds.
    # Ester bond: an oxygen atom linked to a carbonyl carbon.
    ester_smarts = "[OX2][CX3](=O)[#6]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester group(s); expected at least 2 acyl ester bonds"
    
    # 3. Check for diphosphate connectivity.
    # Look for a situation in which two phosphorus atoms (atomic number 15)
    # are bridged by an oxygen.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    diphosphate_found = False
    for atom in p_atoms:
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # oxygen neighbor
                # Check if this oxygen is also attached to a different phosphorus.
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() != atom.GetIdx() and nbr2.GetAtomicNum() == 15:
                        diphosphate_found = True
                        break
            if diphosphate_found:
                break
        if diphosphate_found:
            break
    if not diphosphate_found:
        return False, "Diphosphate bridging (two P atoms connected via an O) not found"
    
    # 4. Check for a diacylglycerol fragment.
    # Genuine CDP-diacylglycerols have a glycerol moiety esterified at two positions.
    # Here we look for a fragment corresponding to:
    #   CH2(OC(=O)[*]) - CH(OC(=O)[*]) - C(OP)
    # which roughly represents a diacylglycerol with the sn-3 position attached to phosphate.
    diacylglycerol_smarts = "[CH2](OC(=O)[*])C(OC(=O)[*])C(OP)"
    diacylglycerol_pattern = Chem.MolFromSmarts(diacylglycerol_smarts)
    if not mol.HasSubstructMatch(diacylglycerol_pattern):
        return False, "Diacylglycerol fragment (glycerol backbone with 2 acyl chains attached and linked to phosphate) not found"
    
    # If all checks pass, classify as CDP-diacylglycerol.
    return True, "Molecule contains a cytidine moiety, diphosphate linkage, acyl ester bonds, and a diacylglycerol fragment consistent with CDP-diacylglycerol"

# Example usage:
if __name__ == "__main__":
    test_smiles = "P(OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O)(OP(OC[C@H](OC(=O)CCCCCCCCCCCCC(C)C)COC(=O)CCCCCCC/C=C\\C=C/CCCCCC)(O)=O)(O)=O"
    result, reason = is_CDP_diacylglycerol(test_smiles)
    print(result, reason)