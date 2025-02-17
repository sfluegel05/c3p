"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CDP-diacylglycerol
Definition: A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) at the 1- and 2-positions.
Improved heuristics:
  - Checks for a cytidine moiety.
  - Checks for a diphosphate linkage (two phosphorus atoms bridged by an oxygen).
  - Checks for at least two acyl ester bonds.
  - Additionally verifies presence of a glycerol backbone typical for diacylglycerol.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol is defined as a CDP-glycerol (i.e. a cytidine diphosphate linked to glycerol)
    with acyl groups (usually fatty acyl chains) at the 1- and 2-positions of the glycerol.
    
    This function applies a four-step heuristic:
      1. Verify the presence of a cytidine moiety.
      2. Check for at least two ester bonds (acyl ester bonds).
      3. Look for a diphosphate connectivity (two phosphorus atoms bridged by an oxygen).
      4. Ensure a glycerol backbone is present.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a CDP-diacylglycerol, False otherwise.
        str: A message describing the classification reason.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Verify the presence of a cytidine moiety.
    # This SMARTS aims to capture a cytosine ring (as found in cytidine).
    cytidine_smarts = "n1ccc(N)nc1=O"
    cytidine_pattern = Chem.MolFromSmarts(cytidine_smarts)
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "Cytidine moiety not found"
    
    # 2. Check for acyl ester bonds.
    # This pattern matches an ester bond: an oxygen atom linked to a carbonyl carbon.
    ester_smarts = "[OX2][CX3](=O)[#6]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester group(s); expected at least 2 acyl ester bonds"

    # 3. Check for diphosphate connectivity.
    # Look for two phosphorus atoms bridged by at least one oxygen.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    diphosphate_found = False
    for atom in p_atoms:
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # oxygen neighbor
                # Check if this oxygen connects to a different phosphorus atom.
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
    
    # 4. Check for a glycerol backbone.
    # A typical glycerol backbone can be heuristically represented by a three-carbon chain bearing hydroxyls.
    # Here we use the SMARTS "OCC(O)CO" to capture a fragment similar to glycerol.
    glycerol_smarts = "OCC(O)CO"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not found"
    
    # If all the above checks pass, we classify the molecule as CDP-diacylglycerol.
    return True, "Molecule contains a cytidine moiety, diphosphate linkage, acyl ester bonds, and a glycerol backbone consistent with CDP-diacylglycerol"

# Example usage:
if __name__ == "__main__":
    test_smiles = "P(OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O)(OP(OC[C@H](OC(=O)CCCCCCCCCCCCC(C)C)COC(=O)CCCCCCC/C=C\\C=C/CCCCCC)(O)=O)(O)=O"
    result, reason = is_CDP_diacylglycerol(test_smiles)
    print(result, reason)