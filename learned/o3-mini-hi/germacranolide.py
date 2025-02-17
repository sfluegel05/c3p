"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: Germacranolide – a sesquiterpene lactone based on the germacrane skeleton.
A germacranolide is defined here by the presence of (1) a lactone ring (cyclic ester) and 
(2) a 10‐membered ring (germacrane skeleton) with mostly carbons and at least two double bonds 
(in order to capture the unsaturated (terpene) nature). Also the overall molecular weight is checked 
to be within a range expected for sesquiterpenes.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide (a sesquiterpene lactone with a germacrane skeleton)
    based on its SMILES string.
    
    The criteria are:
      1. The molecule must be valid.
      2. It must contain a lactone ring: a cyclic ester (the carbonyl carbon and the ester oxygen 
         should both be in a ring).
      3. It must contain a 10-membered ring that has at least 8 carbon atoms and at least 2 non-aromatic double bonds.
         (This extra unsaturation criterion is intended to capture the typical unsaturated germacrane skeleton).
      4. The overall molecular weight should be roughly within the range for sesquiterpene lactones (200 – 400 Da).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a germacranolide, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Check for a lactone ring.
    # Lactone is defined as a cyclic ester, i.e. a carbonyl (C(=O)) linked to an oxygen where both atoms are in a ring.
    lactone_smarts = "[#6;R](=O)[O;R]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring (cyclic ester) found"
    
    # Criterion 2: Search for a 10-membered ring that may represent the germacrane skeleton.
    # We require that: 
    #   a) The ring has 10 atoms.
    #   b) It has at least 8 carbon atoms.
    #   c) It contains at least 2 non-aromatic double bonds (unsaturation typical for terpenes).
    ring_info = mol.GetRingInfo()
    germacrane_found = False
    for ring in ring_info.AtomRings():
        if len(ring) == 10:
            # Count number of carbons.
            carbon_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if carbon_count < 8:
                continue
            # Count double bonds (non-aromatic) among bonds connecting ring atoms.
            double_bond_count = 0
            ring_set = set(ring)
            for bond in mol.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                # Both atoms in the ring?
                if a1 in ring_set and a2 in ring_set:
                    if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
                        double_bond_count += 1
            if double_bond_count >= 2:
                germacrane_found = True
                break
    if not germacrane_found:
        return False, "No 10-membered (germacrane-type) ring with sufficient unsaturation found"
    
    # Criterion 3: Check overall molecular weight as a crude filter.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 400:
        return False, "Molecular weight not in range for typical sesquiterpene lactones"
    
    # If all criteria are met, classify as germacranolide.
    return True, "Contains a lactone ring and a 10-membered unsaturated germacrane skeleton typical for germacranolide"

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES (Provincialin)
    test_smiles = r"C\C=C(/CO)C(=O)OC\C(=C/CO)C(=O)O[C@@H]1C\C(C)=C\C[C@H](OC(C)=O)\C(C)=C/[C@H]2OC(=O)C(=C)[C@H]12"
    result, reason = is_germacranolide(test_smiles)
    print(f"Result: {result}\nReason: {reason}")