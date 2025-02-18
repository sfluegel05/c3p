"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: Nonclassic Icosanoid
Definition: Any biologically active signalling molecule made by oxygenation of C20 fatty acids
            other than the classic icosanoids (the leukotrienes and the prostanoids).

Improved heuristic criteria:
  - The molecule should have roughly 20 carbons (we require between 18 and 23 carbon atoms).
  - Must have at least 4 oxygen atoms.
  - Must have at least 3 double bonds.
  - Must have at least one carbonyl group belonging to either a free carboxylic acid or an ester.
  - Molecules containing benzene rings (aromatic 6‐membered rings) are rejected.
  - Molecules containing nonaromatic 5‐membered rings are rejected.
  - Molecules containing more than one ring are flagged (icosanoid derivatives are usually acyclic
    or contain at most one small heterocycle).
  - Allow up to 3 nitrogen atoms.
  - Overall molecular weight must be below 600 Da.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.

    A nonclassic icosanoid is taken here as a C20 fatty acid derivative that has undergone oxygenation,
    but is not in the classic prostanoid/leukotriene families. The following heuristics are applied:
      - The molecule must have between 18 and 23 carbon atoms (to capture the C20 fatty acid skeleton).
      - The molecule must have at least 4 oxygen atoms.
      - The molecule must have at least 3 double bonds.
      - The molecule must contain at least one carboxyl or ester carbonyl group.
      - Molecules containing benzene rings (fully aromatic 6‐membered rings) are rejected.
      - Molecules containing nonaromatic 5‐membered rings are rejected.
      - Molecules with more than one ring (i.e. multiple cyclic systems) are also rejected.
      - The molecule is allowed up to 3 nitrogen atoms.
      - The overall molecular weight must be less than 600 Da.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a nonclassic icosanoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms (atomic number 6)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 23):
        return False, f"Molecule has {c_count} carbons; expected between 18 and 23 for a C20 fatty acid derivative."

    # Count oxygen atoms (atomic number 8)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Molecule has {o_count} oxygen atoms; expected at least 4 for sufficient oxygenation."
    
    # Count nitrogen atoms (atomic number 7) -- allow up to 3
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count > 3:
        return False, f"Molecule contains {n_count} nitrogen atoms; too many for a fatty acid derivative."

    # Count double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 3:
        return False, f"Found {double_bond_count} double bonds; expected at least 3 for a polyunsaturated fatty acid derivative."
    
    # Check for a carboxyl or ester carbonyl group.
    # SMARTS for carboxylic acid (or deprotonated form) and a common ester pattern.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,O-]")
    ester_pattern = Chem.MolFromSmarts("C(=O)O[!#6]")  # roughly matches an ester if the oxygen is not bonded to another carbon (avoid matching simple acids)
    if not (mol.HasSubstructMatch(acid_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No carboxyl or ester carbonyl group found; expected in an icosanoid‐type molecule."

    # Reject molecules with benzene rings (fully aromatic 6-membered rings)
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if mol.HasSubstructMatch(benzene):
        return False, "Molecule contains benzene rings, which are atypical for icosanoid derivatives."

    # Reject molecules with nonaromatic 5-membered rings.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            aromatic_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
            if aromatic_count == 0:
                return False, "Molecule contains a nonaromatic 5-membered ring (e.g. cyclopentane), indicative of classic prostanoids."

    # Reject molecules with more than one ring (icosanoids are typically acyclic or nearly acyclic)
    if ring_info.NumRings() > 1:
        return False, f"Molecule contains {ring_info.NumRings()} rings; expected at most one ring for an icosanoid derivative."
    
    # Check overall molecular weight (<600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt >= 600:
        return False, f"Molecular weight is {mol_wt:.1f} Da, which is too high for a typical nonclassic icosanoid."
    
    return True, "Molecule is an oxygenated C20 fatty acid derivative (nonclassic icosanoid) lacking classic prostanoid/leukotriene features."

# Example usage (for testing):
if __name__ == "__main__":
    # Test one of the provided examples (20-hydroxylipoxin A4).
    test_smiles = "C(/C=C\\C=C\\[C@H](CCCCCO)O)=C\\C=C\\[C@H]([C@H](CCCC(=O)O)O)O"
    result, reason = is_nonclassic_icosanoid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")