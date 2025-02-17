"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: Germacranolide – a sesquiterpene lactone based on the germacrane skeleton.
A germacranolide is defined here by the presence of:
  (1) a lactone ring (a cyclic ester; i.e. the carbonyl carbon and the ester oxygen share a ring),
  (2) a 10‐membered ring (or an SSSR cycle of 10 atoms) that is predominantly carbon (at least 8 carbons)
      and shows some unsaturation (at least one non‐aromatic double bond) – capturing the terpene nature,
  (3) a molecular weight within the range expected for sesquiterpene lactones (here: 150 – 600 Da),
  (4) and a skeleton containing at least 15 carbon atoms.
Note: This heuristic rule‐based approach is an attempt to balance detection of known examples and rejection
of similar structures that are not true germacranolides.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide (a sesquiterpene lactone with a germacrane skeleton)
    based on its SMILES string.

    The criteria are:
      1. The molecule must be a valid RDKit molecule.
      2. It must contain a lactone ring (a cyclic ester, detected via a SMARTS pattern).
      3. It must contain at least one 10-membered ring (according to the SSSR) with:
            - at least 8 of its atoms being carbons
            - at least one non-aromatic double bond (as a simple proxy for terpene unsaturation)
      4. The overall molecular weight should be within 150 – 600 Da.
      5. The molecule should contain at least 15 carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if classified as a germacranolide, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Check for a lactone ring.
    # Lactone is defined as a cyclic ester: a carbonyl carbon (sp2 C) connected to an oxygen,
    # with both atoms in a ring.
    lactone_smarts = "[#6;R](=O)[O;R]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring (cyclic ester) found"
    
    # Criterion 2: Check for a 10-membered ring with at least 8 carbons and at least one non-aromatic double bond.
    # We use RDKit’s GetSymmSSSR to get a set of rings.
    rings = list(Chem.GetSymmSSSR(mol))
    germacrane_found = False
    for ring in rings:
        if len(ring) == 10:
            # Count carbon atoms in the ring.
            carbon_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if carbon_count < 8:
                continue
            # Count non-aromatic double bonds between atoms in the ring.
            ring_set = set(ring)
            double_bond_count = 0
            # Iterate over bonds; count if both ends are in the ring.
            for bond in mol.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                if a1 in ring_set and a2 in ring_set:
                    if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
                        double_bond_count += 1
            if double_bond_count >= 1:
                germacrane_found = True
                break
    if not germacrane_found:
        return False, "No 10-membered (germacrane-type) ring with sufficient unsaturation found"
    
    # Criterion 3: Check overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 600:
        return False, "Molecular weight not in range for typical sesquiterpene lactones"
    
    # Criterion 4: Check that the molecule has at least 15 carbon atoms.
    carbon_total = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_total < 15:
        return False, "Too few carbons to be a sesquiterpene (expected ≈15 in the core skeleton)"
    
    # If all criteria are met, classify as germacranolide.
    return True, "Contains a lactone ring and a 10-membered unsaturated (germacrane) skeleton typical for germacranolide"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example: Provincialin (should be classified as germacranolide)
    test_smiles = r"C\C=C(/CO)C(=O)OC\C(=C/CO)C(=O)O[C@@H]1C\C(C)=C\C[C@H](OC(C)=O)\C(C)=C/[C@H]2OC(=O)C(=C)[C@H]12"
    result, reason = is_germacranolide(test_smiles)
    print(f"Result: {result}\nReason: {reason}")