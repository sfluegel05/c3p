"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: Nonclassic Icosanoid
Definition: Any biologically active signalling molecule made by oxygenation of C20 fatty acids
            other than the classic icosanoids (the leukotrienes and the prostanoids).

Heuristic criteria used in this classifier (improved):
  - Must have at least 20 carbon atoms.
  - Must have at least 4 oxygen atoms.
  - Must have at least 3 double bonds.
  - Must contain at least one carboxyl AND/OR ester carbonyl group.
    (Both free carboxylic acids “–C(=O)O” and esters “–C(=O)OR” are recognized.)
  - Molecules containing a nonaromatic 5‐membered ring (e.g. a cyclopentane) are rejected.
  - Molecules containing benzene rings (fully aromatic 6‐membered rings) are rejected.
  - Molecules containing more than 2 nitrogen atoms are rejected.
  - The overall molecular weight must be less than 600 Da.
  - Molecules with multiple (more than one) rings are also flagged, as icosanoid derivatives are usually
    acyclic or contain at most one small heterocycle.
    
Note: These rules are heuristic; more refinements may be needed for some borderline cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.

    A nonclassic icosanoid is defined here as an oxygenated C20 fatty acid derivative that does not 
    belong to the classic prostanoid/leukotriene family. The following heuristics are applied:
      - Must have at least 20 carbon atoms.
      - Must have at least 4 oxygen atoms.
      - Must have at least 3 double bonds.
      - Must contain at least one carboxyl (or ester) carbonyl group.
      - Must not contain a nonaromatic 5‐membered ring (which is typical of prostanoids).
      - Must not contain any benzene rings.
      - Allows a maximum of 2 nitrogen atoms.
      - Should not contain multiple rings (more than one).
      - Overall molecular weight must be less than 600 Da.
      
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the compound is classified as a nonclassic icosanoid, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms (atomic number 6)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Molecule has {c_count} carbon atoms; expected at least 20 for a C20 fatty acid derivative."
    
    # Count oxygen atoms (atomic number 8)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Molecule has {o_count} oxygen atoms; expected at least 4 for sufficient oxygenation."
    
    # Count nitrogen atoms (atomic number 7), allow at most 2
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count > 2:
        return False, f"Molecule contains {n_count} nitrogen atoms; too many for a fatty acid derivative."
    
    # Count double bonds (unsaturation)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 3:
        return False, f"Found {double_bond_count} double bonds; expected at least 3 for a polyunsaturated fatty acid derivative."
    
    # Check for the presence of a carboxyl or ester carbonyl group.
    # Define SMARTS: one for carboxylic acid (including deprotonated form) and one for esters.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,O-]")  # matches -C(=O)OH or its deprotonated form
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")       # rough pattern for esters (will match -C(=O)OR)
    if not (mol.HasSubstructMatch(acid_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No carboxyl or ester carbonyl group found; expected in icosanoid‐type molecules."
    
    # Exclude molecules with a benzene ring (fully aromatic 6-membered ring)
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if mol.HasSubstructMatch(benzene):
        return False, "Molecule contains benzene rings, which are atypical for icosanoid derivatives."
    
    # Exclude molecules containing a nonaromatic 5-membered ring.
    # Loop over the rings in the molecule.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            # Count aromatic atoms in this ring.
            aromatic_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
            # If none of the atoms are aromatic, then reject this ring.
            if aromatic_count == 0:
                return False, "Molecule contains a nonaromatic 5-membered ring (cyclopentane), indicative of classic prostanoids."
    
    # Optionally, reject molecules with multiple rings (more than one)
    if ring_info.NumRings() > 1:
        return False, f"Molecule contains {ring_info.NumRings()} rings; expected at most one ring for a typical icosanoid derivative."
    
    # Check overall molecular weight (<600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 600:
        return False, f"Molecular weight is {mol_wt:.1f} Da, which is too high for a typical nonclassic icosanoid."
    
    # Passed all heuristic checks
    return True, "Molecule is an oxygenated C20 fatty acid derivative (nonclassic icosanoid) lacking classic prostanoid/leukotriene features."

# Example usage (for testing):
if __name__ == "__main__":
    # Test with one example from the provided list (20-hydroxylipoxin A4).
    test_smiles = "C(/C=C\\C=C\\[C@H](CCCCCO)O)=C\\C=C\\[C@H]([C@H](CCCC(=O)O)O)O"
    result, reason = is_nonclassic_icosanoid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")