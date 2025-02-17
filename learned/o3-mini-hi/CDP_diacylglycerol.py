"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
#!/usr/bin/env python
"""
Classifies: CDP-diacylglycerol
Definition: A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) at the 1- and 2-positions.
Improved heuristics:
  - Check for a cytidine moiety using a SMARTS pattern.
  - Require a strict diacylglycerol fragment—with two acyl ester groups on a three-carbon glycerol backbone 
    (i.e. CH2(OC(=O)[#6])-[CH](OC(=O)[#6])-[CH2]OP).
  - Verify diphosphate connectivity, meaning two phosphorus atoms connected via a bridging oxygen.
"""

from rdkit import Chem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    
    The heuristic criteria are as follows:
      1. The molecule must contain a cytidine moiety, captured by the SMARTS "n1ccc(N)nc1=O".
      2. The molecule must have a diacylglycerol fragment that looks like:
         CH2(OC(=O)[#6])-[CH](OC(=O)[#6])-[CH2]OP
         which represents a glycerol backbone where positions 1 and 2 are acylated (esterified) 
         and the third carbon is linked (via an O) to a phosphate.
      3. The molecule must exhibit diphosphate connectivity; that is, two phosphorus atoms (atomic number 15)
         that are bridged by an oxygen.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a CDP-diacylglycerol, False otherwise.
        str: A message explaining the reasoning behind the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for cytidine moiety.
    cytidine_smarts = "n1ccc(N)nc1=O"
    cytidine_pattern = Chem.MolFromSmarts(cytidine_smarts)
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "Cytidine moiety not found"
    
    # 2. Check for a strict diacylglycerol fragment.
    # This pattern demands a 3-carbon (glycerol) backbone with esterified acyl groups
    # at the first two carbons and a phosphate connected at the third carbon.
    diacylglycerol_smarts = "[CH2](OC(=O)[#6])-[CH](OC(=O)[#6])-[CH2]OP"
    diacylglycerol_pattern = Chem.MolFromSmarts(diacylglycerol_smarts)
    if not mol.HasSubstructMatch(diacylglycerol_pattern):
        return False, "Diacylglycerol fragment (glycerol backbone with two acyl esters and phosphate linkage) not found"
    
    # 3. Check for diphosphate connectivity:
    # Two phosphorus atoms connected via a bridging oxygen.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    diphosphate_found = False
    for atom in p_atoms:
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # oxygen neighbor
                # Check whether this oxygen also connects to a different phosphorus atom.
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() != atom.GetIdx() and nbr2.GetAtomicNum() == 15:
                        diphosphate_found = True
                        break
            if diphosphate_found:
                break
        if diphosphate_found:
            break
    if not diphosphate_found:
        return False, "Diphosphate bridging (two phosphorus atoms connected via an oxygen) not found"
    
    return True, "Molecule contains a cytidine moiety, diphosphate linkage, and a diacylglycerol fragment with two acyl ester bonds consistent with CDP-diacylglycerol"

# Example usage:
if __name__ == "__main__":
    test_smiles = "P(OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O)(OP(OC[C@H](OC(=O)CCCCCCCCCCCCC(C)C)COC(=O)CCCCCCC/C=C\\C=C/CCCCCC)(O)=O)(O)=O"
    result, reason = is_CDP_diacylglycerol(test_smiles)
    print(result, reason)