"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: Nonclassic Icosanoid
Definition: Any biologically active signalling molecule made by oxygenation of C20 fatty acids
            other than the classic icosanoids (the leukotrienes and the prostanoids).

Heuristic criteria used in this classifier:
  - The molecule must have at least 20 carbon atoms.
  - It should have at least 4 oxygen atoms.
  - It should contain at least 3 double bonds (reflecting polyunsaturation).
  - There should be at least one carboxyl or ester carbonyl group (indicative of a fatty acid end).
  - Molecules with a non‐aromatic 5‐membered ring (e.g. a cyclopentane ring) are rejected.
  - Molecules with benzene rings (fully aromatic 6‐membered rings) are rejected.
  - Molecules containing more than 2 nitrogen atoms are rejected.
  - The overall molecular weight should be less than 600 Da.
  
Note: These rules are heuristic and based on our analysis of true and false positives;
      they may need further refinement.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.

    A nonclassic icosanoid is defined as an oxygenated C20 fatty acid derivative that does 
    not belong to the classic (leukotriene/prostanoid) family. The following heuristic are applied:
      - Must have at least 20 carbon atoms.
      - Must have at least 4 oxygen atoms.
      - Must have at least 3 double bonds.
      - Must contain at least one carboxyl (or ester) carbonyl group.
      - Must not contain a nonaromatic 5‐membered ring (cyclopentane ring typical of prostanoids).
      - Must not contain any benzene rings.
      - Allows at most 2 nitrogen atoms.
      - Overall mol. weight must be less than 600 Da.
      
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the compound is classified as a nonclassic icosanoid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Molecule has {c_count} carbons; expected at least 20 for a C20 fatty acid derivative."
    
    # Count oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Molecule has only {o_count} oxygens; insufficient oxygenation for an icosanoid derivative."
    
    # Count nitrogen atoms. Allow up to 2.
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count > 2:
        return False, f"Molecule contains {n_count} nitrogen atoms; too many for a fatty acid derivative."
    
    # Count double bonds (unsaturation)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 3:
        return False, f"Only {double_bond_count} double bonds found; expected at least 3 for a polyunsaturated fatty acid."
    
    # Check for presence of a carboxyl or ester carbonyl group.
    # SMARTS for a carbonyl with an -OH or -OR: [CX3](=O)[OX1H0,OX2H1]
    carbonyl_smarts = Chem.MolFromSmarts("[CX3](=O)[OX1,H]")
    if not mol.HasSubstructMatch(carbonyl_smarts):
        return False, "No carboxyl/ester carbonyl group found; expected in icosanoid‐type molecules."
    
    # Exclude molecules containing a nonaromatic 5-membered ring (e.g. cyclopentane)
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == 5:
            # Count aromatic atoms in the ring.
            aromatic_atoms = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
            if aromatic_atoms == 0:
                return False, "Molecule contains a nonaromatic 5-membered ring (cyclopentane), indicative of prostanoids."
    
    # Exclude molecules that contain benzene rings 
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if mol.HasSubstructMatch(benzene):
        return False, "Molecule contains benzene rings, which are atypical for icosanoid derivatives."
    
    # Check overall molecular weight (icosanoids are usually not very heavy)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 600:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too high for a typical nonclassic icosanoid."
    
    # Passed all heuristic checks
    return True, "Molecule is an oxygenated C20 fatty acid derivative (nonclassic icosanoid) lacking classic prostanoid/leukotriene features."

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = "C(/C=C\\C=C\\[C@H](CCCCCO)O)=C\\C=C\\[C@H]([C@H](CCCC(=O)O)O)O"  # Example: 20-hydroxylipoxin A4
    result, reason = is_nonclassic_icosanoid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")