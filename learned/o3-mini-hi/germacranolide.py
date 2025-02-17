"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: Germacranolide – a sesquiterpene lactone based on the germacrane skeleton.
A germacranolide is defined here by the following criteria:
  (1) It contains a lactone ring (a cyclic ester where the carbonyl carbon and the ester oxygen are in a ring).
  (2) It has at least one ring of size between 10 and 12 atoms that is predominately carbon (>=80% of its atoms are C)
      and has at least one non‐aromatic double bond within that ring.
      (This aims to capture the germacrane skeleton even when it appears as a fused ring system.)
  (3) Its molecular weight is within the range typical for sesquiterpene lactones (150–600 Da).
  (4) The overall molecule contains at least 15 carbon atoms.
Note: This heuristic approach is not perfect; some true germacranolides may be missed and some non‐germacranolides
may still pass. The goal is to improve the F1 score relative to the previous attempt.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide (a sesquiterpene lactone with a germacrane skeleton)
    based on its SMILES string.
    
    The improved criteria are:
      1. The molecule must be parsed correctly.
      2. It must contain a lactone ring (cyclic ester).
      3. It must contain at least one ring (from those reported in GetRingInfo().AtomRings())
         whose size is between 10 and 12 atoms, where at least 80% of the atoms are carbons and
         which has at least one non‐aromatic double bond within that ring.
      4. The overall molecular weight should be within 150 – 600 Da.
      5. The molecule must have at least 15 carbon atoms.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a germacranolide, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Check for a lactone ring.
    # A lactone is defined here as a cyclic ester: a carbonyl carbon (sp2 C) connected to an oxygen all within a ring.
    lactone_smarts = "[#6;R](=O)[OX2;R]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring (cyclic ester) found"
    
    # Criterion 2: Look for a candidate germacrane-type ring.
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_ring_found = False
    for ring in ring_info:
        ring_size = len(ring)
        # Accept rings that have size between 10 and 12 atoms:
        if not (10 <= ring_size <= 12):
            continue
        # Count how many atoms in the ring are carbons.
        carbon_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        # Must have at least 80% carbons:
        if carbon_count / ring_size < 0.8:
            continue
        
        # Count non-aromatic double bonds within this ring.
        # (Only count bonds where both endpoints are in the ring.)
        ring_set = set(ring)
        double_bond_count = 0
        # Iterate over bonds in the molecule and count those within the ring having a double bond.
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_set and a2 in ring_set:
                # Check for double bond and ensure it is not aromatic.
                if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
                    double_bond_count += 1
        if double_bond_count >= 1:
            candidate_ring_found = True
            break  # Accept the first candidate ring that fulfills the condition.
    
    if not candidate_ring_found:
        return False, "No 10–12-membered ring with sufficient carbon content and unsaturation found"
    
    # Criterion 3: Check overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 600:
        return False, "Molecular weight not in range for typical sesquiterpene lactones (150–600 Da)"
    
    # Criterion 4: Check that the molecule has at least 15 carbon atoms.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 15:
        return False, "Too few carbons (<15) for a sesquiterpene lactone skeleton"
    
    return True, "Contains a lactone ring and a 10–12-membered, predominantly carbon, unsaturated ring typical for germacranolide"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example: Provincialin (should be classified as germacranolide)
    test_smiles = r"C\C=C(/CO)C(=O)OC\C(=C/CO)C(=O)O[C@@H]1C\C(C)=C\C[C@H](OC(C)=O)\C(C)=C/[C@H]2OC(=O)C(=C)[C@H]12"
    result, reason = is_germacranolide(test_smiles)
    print(f"Input SMILES: {test_smiles}")
    print(f"Result: {result}\nReason: {reason}")