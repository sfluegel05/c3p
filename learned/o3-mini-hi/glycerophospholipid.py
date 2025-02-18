"""
Classifies: CHEBI:37739 glycerophospholipid
"""
#!/usr/bin/env python
"""
Classifies: glycerophospholipid 
Definition:
  A glycerophospholipid is any glycerolipid having a phosphate group ester‐linked to 
  a terminal (primary: CH2) carbon of a glycerol backbone plus at least one acyl ester 
  (OC(=O)) group.
Heuristic criteria in this revision:
  1. The molecule must be ≥400 Da (to avoid many small phosphate‐containing compounds).
  2. It must contain a glycerol backbone. We search for a three‐carbon contiguous chain 
     with two terminal CH2 groups. Then at least one terminal CH2 must carry an oxygen 
     which is directly bonded to a phosphorus (i.e. an ester linkage to phosphate).
  3. There must be at least one acyl ester group (OC(=O)) elsewhere in the molecule 
     (i.e. not counting oxygens directly attached to phosphorus).
Because of the chemical diversity, the filter is heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    Heuristically, this function checks that:
      (a) the molecule has a molecular weight ≥400 Da,
      (b) it contains a glycerol backbone (a three‐carbon chain with terminal CH2 groups)
          and one of the terminal CH2 carbons is esterified to a phosphate (i.e. O–P bonding),
      (c) at least one acyl ester group (OC(=O)) is present outside the phosphate motif.
    
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule is classified as a glycerophospholipid; False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Check molecular weight ≥ 400 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a glycerophospholipid"
    
    # Criterion 2: Look for a glycerol backbone.
    # We search for a three-carbon chain with pattern: terminal CH2 - CH - terminal CH2.
    glycerol_smarts = "[CH2][CH][CH2]"
    glycerol_pat = Chem.MolFromSmarts(glycerol_smarts)
    backbone_matches = mol.GetSubstructMatches(glycerol_pat)
    if not backbone_matches:
        return False, "No glycerol backbone (three contiguous carbon chain with terminal CH2 groups) found"
    
    # For each glycerol backbone match, check if one of the terminal (first or third) atoms 
    # is connected via an oxygen to a phosphorus.
    phosphate_attached = False
    for match in backbone_matches:
        # match is a tuple of three atom indices in order.
        for pos in (0, 2):  # check each terminal CH2
            carbon = mol.GetAtomWithIdx(match[pos])
            # Check that this carbon really is CH2: exactly 2 hydrogens.
            if carbon.GetTotalNumHs() != 2:
                continue
            # Look at neighbors of this carbon that are oxygen.
            for nbr in carbon.GetNeighbors():
                if nbr.GetSymbol() != 'O':
                    continue
                # To be a phosphate ester linkage, that oxygen must be bonded to a P.
                for o_nbr in nbr.GetNeighbors():
                    if o_nbr.GetSymbol() == 'P':
                        phosphate_attached = True
                        break
                if phosphate_attached:
                    break
            if phosphate_attached:
                break
        if phosphate_attached:
            break
    if not phosphate_attached:
        return False, "No terminal CH2 of a glycerol backbone found that is ester-linked via oxygen to phosphorus"
    
    # Criterion 3: Look for at least one acyl ester group (OC(=O)) that is not directly attached to a phosphorus.
    # The SMARTS "[O;!$(*~P)]C(=O)" tries to match an oxygen (not attached to any phosphorus)
    # singly bonded to a carbonyl carbon.
    ester_smarts = "[O;!$(*~P)]C(=O)"
    ester_pat = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pat)
    if len(ester_matches) < 1:
        return False, "No acyl ester (OC(=O)) group found that would indicate a fatty acyl chain"
    
    return True, ("Molecule has a molecular weight ≥400 Da, a glycerol backbone with a terminal CH2 "
                  "ester-linked via oxygen to phosphate, and at least one acyl ester group, consistent with a "
                  "glycerophospholipid structure.")

# Example usage (for testing)
if __name__ == "__main__":
    # A known glycerophospholipid SMILES example:
    test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)(O)(O)=O"  # 1-(11Z,14Z-eicosadienoyl)-glycero-3-phosphate
    result, reason = is_glycerophospholipid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)