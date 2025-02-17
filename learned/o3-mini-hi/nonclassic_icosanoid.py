"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: Nonclassic Icosanoid
Definition: Any biologically active signalling molecule made by oxygenation of C20 fatty acids
            other than the classic icosanoids (the leukotrienes and the prostanoids).

Improved heuristic criteria for classification (heuristic only):
  - Total number of carbon atoms must be between 19 and 32.
  - The molecule must contain at least 4 oxygen atoms.
  - The oxygen:carbon ratio must be ≥ 0.20.
  - The molecule must have at least 3 double bonds.
  - The molecule should display at least one conjugated diene/triene motif 
    (detected via a SMARTS for consecutive double bonds).
  - The molecule must have at least one carbonyl group as part of either a carboxylic acid
    (SMARTS: "C(=O)[O;H,-]") or an ester (SMARTS: "C(=O)OC", which now catches methyl esters).
  - Molecules containing benzene rings (fully aromatic 6‐membered rings) are rejected.
  - Molecules that contain nonaromatic 5‐membered rings are rejected.
  - Molecules with more than one ring (i.e. multiple cyclic systems) are rejected.
  - Up to 3 nitrogen atoms are allowed.
  - Overall molecular weight must be below 600 Da.
  
Note: This is a heuristic approach – subtle cases (and many edge cases) may be mis‐classified.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.

    The heuristics applied are:
      - Total carbon count between 19 and 32 (allowing modifications to a C20 fatty acid).
      - At least 4 oxygen atoms and an O/C ratio >= 0.20.
      - At least 3 double bonds.
      - At least one instance of conjugated double bonds (i.e. a pattern like C=C-C=C).
      - At least one carbonyl group as part of either a carboxylic acid or an ester.
      - Molecule must not contain benzene rings.
      - Molecule must not contain nonaromatic cyclopentane (5-membered) rings.
      - At most one ring overall.
      - At most 3 nitrogen atoms.
      - Molecular weight below 600 Da.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a nonclassic icosanoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Count carbon, oxygen, and nitrogen atoms
    atoms = list(mol.GetAtoms())
    c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 7)
    
    # Check carbon count (expected to derive from a C20 fatty acid, allow for modifications)
    if not (19 <= c_count <= 32):
        return False, f"Molecule has {c_count} carbons; expected between 19 and 32 for a C20-derived molecule."
    
    # Check oxygen count (need sufficient oxygenation) and oxygen/carbon ratio
    if o_count < 4:
        return False, f"Molecule has {o_count} oxygen atoms; expected at least 4 for oxygenation."
    if (o_count / c_count) < 0.20:
        return False, f"Oxygen/carbon ratio ({o_count}/{c_count} ≈ {(o_count/c_count):.2f}) is too low."
    
    # Allow at most 3 nitrogen atoms
    if n_count > 3:
        return False, f"Molecule has {n_count} nitrogen atoms; too many for a fatty acid derivative."
    
    # Count double bonds (only count bonds that are double bonds and between heavy atoms)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 3:
        return False, f"Found {double_bond_count} double bonds; expected at least 3 for a polyunsaturated derivative."
    
    # Check for conjugated double bonds: search for a pattern [C]=[C]-[C]=[C]
    conjugated_pattern = Chem.MolFromSmarts("[#6]=[#6]-[#6]=[#6]")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "No conjugated double bond pattern (e.g. C=C-C=C) found; expected in an icosanoid-type molecule."
    
    # Check for carbonyl presence as part of a carboxylic acid or an ester.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,O-]")
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")  # now matches a general ester (e.g. methyl esters)
    if not (mol.HasSubstructMatch(acid_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No carboxyl or ester carbonyl group found; expected in an icosanoid-type molecule."
    
    # Reject molecules with benzene rings (aromatic, six-membered rings)
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if mol.HasSubstructMatch(benzene):
        return False, "Molecule contains benzene rings, atypical for an icosanoid derivative."
    
    # Use ring information to reject nonaromatic 5-membered rings and molecules with >1 ring.
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:
        return False, f"Molecule contains {ring_info.NumRings()} rings; expected at most one ring for an icosanoid derivative."
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            # If none of the atoms in the ring is aromatic then reject (nonaromatic 5-membered ring)
            if not any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                return False, "Molecule contains a nonaromatic 5-membered ring, which is atypical for nonclassic icosanoids."
    
    # Check overall molecular weight (<600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt >= 600:
        return False, f"Molecular weight is {mol_wt:.1f} Da, which is too high for a typical nonclassic icosanoid."
    
    return True, "Molecule is an oxygenated C20 fatty acid derivative (nonclassic icosanoid) lacking classic prostanoid/leukotriene features."


# Example usage (for testing purposes):
if __name__ == "__main__":
    # Use one of the reported examples (20-hydroxylipoxin A4)
    test_smiles = "C(/C=C\\C=C\\[C@H](CCCCCO)O)=C\\C=C\\[C@H]([C@H](CCCC(=O)O)O)O"
    result, reason = is_nonclassic_icosanoid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")