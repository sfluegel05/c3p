"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-5β-steroid
Definition: Any 3-oxo steroid that has beta- configuration at position 5.
This improved implementation uses heuristics:
  • It requires that the molecule possesses a fused, non‐aromatic ring system 
    containing at least 17 carbon atoms (the steroid nucleus).
  • It requires a ring-bound ketone group (the 3-oxo moiety).
  • It checks that the isomeric SMILES string includes at least one '@@'
    as a proxy for a 5β configuration.
Note: A full analysis would require explicit detection of the steroid numbering system.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5β-steroid.
    
    Heuristic criteria:
      1. The molecule must contain a fused, non-aromatic ring system forming a steroid nucleus.
         Here we require that the union of atoms that participate in non-aromatic rings 
         contains at least 17 carbons.
      2. The molecule must contain at least one ketone group on a ring 
         (the 3-oxo moiety), detected using the SMARTS "[R]C(=O)[R]".
      3. At least one chiral center in the structure is indicated with a '@@' in the isomeric
         SMILES string – used here as a proxy for 5β configuration.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is likely a 3-oxo-5β-steroid, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure that stereochemistry is assigned.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # --- Criterion 1: Steroid nucleus ---
    # We consider only non-aromatic rings (the steroid nucleus is saturated).
    ring_info = mol.GetRingInfo()
    non_aromatic_ring_atoms = set()
    for ring in ring_info.AtomRings():
        # Check that the entire ring is non-aromatic.
        if all(not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            non_aromatic_ring_atoms.update(ring)
    # Count the number of carbon atoms in the collected rings.
    ring_carbon_count = 0
    for idx in non_aromatic_ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            ring_carbon_count += 1
    if ring_carbon_count < 17:
        return False, f"Non-aromatic ring system contains only {ring_carbon_count} carbons; expected at least 17 for a steroid nucleus"
    
    # --- Criterion 2: 3-oxo (ring-bound ketone) ---
    # SMARTS: "[R]C(=O)[R]" matches a ketone with both neighbors being ring members.
    ketone_pattern = Chem.MolFromSmarts("[R]C(=O)[R]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ring-bound ketone (3-oxo) group found"
    
    # --- Criterion 3: Proxy for 5β configuration ---
    # We require that the isomeric SMILES includes at least one chiral specification "@@".
    iso_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
    if "@@" not in iso_smi:
        return False, "No chiral center with '@@' (indicative of beta configuration) detected"
    
    return True, "Molecule has a fused non-aromatic ring system with at least 17 carbons, a ring-bound ketone (3-oxo), and a '@@' chiral center indicative of 5β configuration"

# Example usage (uncomment below lines to test):
#if __name__ == "__main__":
#    test_smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2C(=O)CO"  # 5β-dihydrodeoxycorticosterone
#    result, reason = is_3_oxo_5beta_steroid(test_smiles)
#    print(result, reason)