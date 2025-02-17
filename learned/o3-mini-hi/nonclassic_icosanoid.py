"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: Nonclassic Icosanoid
Definition: Any biologically active signalling molecule made by oxygenation of C20 fatty acids
            other than the classic icosanoids (the leukotrienes and the prostanoids).

Heuristic criteria used in this classifier:
  - The molecule must have at least 20 carbon atoms, indicating a fatty acid backbone.
  - The molecule must be oxygenated – we require at least 4 oxygen atoms so as to capture oxygenation
    beyond the carboxyl group.
  - There must be at least 3 double bonds, reflecting the polyunsaturation typically seen in fatty acids.
  - To avoid classic prostanoids, any non‐aromatic 5‐membered ring is rejected.
  - The molecule must not contain any nitrogen atoms (to avoid peptides and other non‐lipid compounds).
  - We also check for the presence of a carboxyl/ester carbonyl group using a SMARTS pattern.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    
    A nonclassic icosanoid is defined as an oxygenated derivative of a C20 fatty acid,
    excluding the classic icosanoids (leukotrienes and prostanoids).
    
    Heuristics applied:
      - Must have at least 20 carbon atoms.
      - Must have at least 4 oxygen atoms.
      - Must have at least 3 double bonds.
      - Must contain at least one carboxyl or ester carbonyl group.
      - Must not contain any nitrogen atoms.
      - Must not have a non‐aromatic 5‐membered ring (cyclopentane), which is typical of prostanoids.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a nonclassic icosanoid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Molecule has {c_count} carbons; expected at least 20 for a C20 fatty acid derivative."
    
    # Exclude molecules that contain any nitrogen atoms (e.g. peptides).
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count > 0:
        return False, f"Molecule contains {n_count} nitrogen atoms; likely not a lipid-derived icosanoid."
        
    # Ensure sufficient oxygenation.
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Molecule has only {o_count} oxygen atoms; insufficient oxygenation required for an icosanoid."
    
    # Count double bonds (unsaturations).
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 3:
        return False, f"Only {double_bond_count} double bonds found; expected at least 3 for polyunsaturation."
    
    # Check for the presence of a carboxyl or ester carbonyl group.
    # This simple SMARTS looks for C(=O)O (acid or ester) functionality.
    carbonyl_smarts = Chem.MolFromSmarts("[CX3](=O)[OX1,H]")
    if not mol.HasSubstructMatch(carbonyl_smarts):
        return False, "No carboxyl/ester carbonyl group found; expected in icosanoid-type molecules."
    
    # Exclude molecules containing a non-aromatic 5-membered ring 
    # (a cyclopentane ring is characteristic of prostanoids).
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == 5:
            # Count aromatic atoms in the ring.
            aromatic_atoms = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
            if aromatic_atoms == 0:
                return False, "Molecule contains a non-aromatic 5-membered (cyclopentane) ring, indicative of prostanoids."
    
    # If all heuristic tests pass, we classify the molecule as a nonclassic icosanoid.
    return True, "Molecule is an oxygenated fatty acid derivative (nonclassic icosanoid) lacking classic prostanoid/leukotriene features."

# Example usage (you may test with provided SMILES strings):
if __name__ == "__main__":
    test_smiles = "C(/C=C\\C=C\\[C@H](CCCCCO)O)=C\\C=C\\[C@H]([C@H](CCCC(=O)O)O)O"  # 20-hydroxylipoxin A4
    result, reason = is_nonclassic_icosanoid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")