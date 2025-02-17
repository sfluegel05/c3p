"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-5β-steroid
Definition: Any 3-oxo steroid that has beta- configuration at position 5.
NOTE: This implementation uses heuristics:
  • It requires that the molecule contains at least four fused rings (the steroid nucleus).
  • It must contain a ketone group within a ring (for the 3-oxo).
  • It checks that the isomeric SMILES contains at least one chiral center 
    annotated with '@@' (our proxy for beta configuration at C5).
A full assignment of steroid numbering and configuration would require a more elaborate analysis.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5β-steroid.
    
    Heuristic criteria:
      1. The molecule must have at least 4 rings (steroid nucleus).
      2. The molecule must contain at least one ketone group on a ring 
         (i.e. a 3-oxo moiety).
      3. At least one chiral center in the structure (in the fused ring system)
         is indicated with a '@@' in the isomeric SMILES string – used here as a proxy
         for a β- configuration at C5.
    
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
    
    # Ensure that stereochemistry is assigned
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Check for the steroid nucleus: assume at least 4 rings (the steroid tetracyclic system)
    ri = mol.GetRingInfo()
    num_rings = ri.NumRings()
    if num_rings < 4:
        return False, f"Only {num_rings} rings detected; a steroid nucleus requires at least 4 fused rings"
    
    # Look for a ketone group in a ring (3-oxo). We require the carbonyl carbon to be in a ring.
    # The SMARTS pattern "[R]C(=O)[R]" searches for a carbonyl group with both neighbor atoms coming from a ring.
    ketone_pattern = Chem.MolFromSmarts("[R]C(=O)[R]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group (3-oxo) found in a ring"
    
    # As a proxy for the required beta configuration at position 5, 
    # we check that the isomeric SMILES string includes at least one "@@" chiral specification.
    # In many steroid SMILES, a beta-oriented substituent may be indicated using the '@@' notation.
    iso_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
    if "@@" not in iso_smi:
        return False, "No chiral center with '@@' (indicative of beta configuration) detected"
    
    # If all the above criteria are met then we assume the molecule is a 3-oxo-5β-steroid.
    return True, "Molecule has at least 4 fused rings, a ring-bound ketone (3-oxo), and a '@@' chiral center indicative of 5β configuration"

# Example usage (uncomment to test):
#if __name__ == "__main__":
#    test_smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2C(=O)CO"  # 5β-dihydrodeoxycorticosterone
#    result, reason = is_3_oxo_5beta_steroid(test_smiles)
#    print(result, reason)