"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
#!/usr/bin/env python
"""
Classifies: CDP-diacylglycerol
Definition: A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) at the 1- and 2-positions.
Improved heuristics:
  - Check for a cytidine moiety using a SMARTS pattern.
  - Require a diacylglycerol fragment defined as a three-carbon glycerol backbone with 
    esterified acyl groups on the first two carbons and a phosphate (OP) on the third carbon.
    In order to avoid matching cyclic or otherwise misassembled fragments, we require the three 
    backbone carbons to be acyclic (R0).
  - Verify diphosphate connectivity by looking for a single bridging oxygen that connects exactly two phosphorus atoms.
"""

from rdkit import Chem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    
    Criteria:
      1. Must contain a cytidine moiety (SMARTS: "n1ccc(N)nc1=O").
      2. Must contain a diacylglycerol fragment defined as:
           [CH2;R0](OC(=O)[*])-[CH;R0](OC(=O)[*])-[CH2;R0]OP
         which is a glycerol backbone with the correct acylations and phosphate attachment.
      3. Must have diphosphate connectivity, that is, one oxygen atom bridging exactly two phosphorus atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule fits the criteria for CDP-diacylglycerol, False otherwise.
        str: Explanation of the classification decision.
    """
    # Convert SMILES to molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for cytidine moiety.
    cytidine_smarts = "n1ccc(N)nc1=O"
    cytidine_pattern = Chem.MolFromSmarts(cytidine_smarts)
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "Cytidine moiety not found"
    
    # 2. Check for a diacylglycerol fragment.
    # The SMARTS pattern insists on an acyclic glycerol backbone with two esterified acyl groups and phosphate.
    diacylglycerol_smarts = "[CH2;R0](OC(=O)[*])-[CH;R0](OC(=O)[*])-[CH2;R0]OP"
    diacylglycerol_pattern = Chem.MolFromSmarts(diacylglycerol_smarts)
    if not mol.HasSubstructMatch(diacylglycerol_pattern):
        return False, "Diacylglycerol fragment (acyclic glycerol with two acyl esters and phosphate linkage) not found"
    
    # 3. Check for diphosphate connectivity.
    # We search for an oxygen atom with exactly two phosphorus neighbors.
    bridging_oxygen_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # oxygen atom
            # count phosphorus neighbors for this oxygen
            p_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 15]
            if len(p_neighbors) == 2:
                bridging_oxygen_count += 1
    if bridging_oxygen_count != 1:
        return False, "Diphosphate bridging oxygen (one oxygen bonded to exactly two phosphorus atoms) not found"
    
    return True, "Molecule contains a cytidine moiety, proper diacylglycerol fragment, and diphosphate connectivity consistent with CDP-diacylglycerol"

# Example usage:
if __name__ == "__main__":
    test_smiles = "P(OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O)(OP(OC[C@H](OC(=O)CCCCCCCCCCCCC(C)C)COC(=O)CCCCCCC/C=C\\C=C/CCCCCC)(O)=O)(O)=O"
    result, reason = is_CDP_diacylglycerol(test_smiles)
    print(result, reason)