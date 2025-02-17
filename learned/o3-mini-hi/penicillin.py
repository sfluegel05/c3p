"""
Classifies: CHEBI:17334 penicillin
"""
#!/usr/bin/env python
"""
Classifies: Substituted penams (penicillins) defined as:
  “Any member of the group of substituted penams containing two methyl substituents at position 2, 
   a carboxylate substituent at position 3 and a carboxamido group at position 6.”
    
This implementation uses heuristic SMARTS patterns:
  - A simplified penam core pattern capturing the bicyclic β‐lactam/thiazolidine framework.
  - A pattern for a sulfur-bound carbon with two methyl groups.
  - A generic carboxylate fragment.
  - A generic carboxamido fragment.
    
Note: Because of the complexity and variability in penicillin structures, this method is heuristic.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_penicillin(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a penicillin based on:
      - Presence of a fused bicyclic penam core (β‐lactam fused to thiazolidine)
      - Two methyl substituents on the carbon adjacent to sulfur (position 2)
      - A carboxylate substituent (e.g. C(=O)O) presumed at position 3
      - A carboxamido group (e.g. N-C(=O)) presumed at position 6

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a penicillin, False otherwise.
        str: Reason for the classification decision.
    """

    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Heuristic: Look for the penam core.
    # This simplified SMARTS attempts to capture a fused bicyclic system comprising:
    #  - A beta-lactam ring (N1C(=O)...C1=O)
    #  - A thiazolidine ring containing an S atom and a dimethyl-substituted carbon.
    # Note: Stereochemistry is removed.
    penam_core_smarts = "N1C(=O)C2SC(C)(C)C(N2C1=O)"
    penam_core = Chem.MolFromSmarts(penam_core_smarts)
    if not mol.HasSubstructMatch(penam_core):
        return False, "Molecule does not contain the expected penam core"

    # Check for two methyl substituents at position 2.
    # We look for a sulfur (S) attached to a carbon that bears two methyl groups.
    # The pattern [SX2]-C([CH3])([CH3]) will search for an S bound (via a single bond) 
    # to a carbon with two CH3 groups.
    dimethyl_pattern = Chem.MolFromSmarts("[$(S)-C([CH3])([CH3])]")
    if not mol.HasSubstructMatch(dimethyl_pattern):
        return False, "Missing two methyl substituents at position 2"

    # Check for a carboxylate substituent.
    # We use a generic pattern for a carboxylate group:
    #   C(=O)O  (which may occur as an acid or as a deprotonated carboxylate).
    # Note: This pattern also matches esters but is used here as a heuristic.
    carboxylate1 = Chem.MolFromSmarts("C(=O)[O-]")
    carboxylate2 = Chem.MolFromSmarts("C(=O)O")
    if not (mol.HasSubstructMatch(carboxylate1) or mol.HasSubstructMatch(carboxylate2)):
        return False, "Missing a carboxylate substituent (C(=O)O group)"

    # Check for a carboxamido group.
    # We search for an N-C(=O) fragment.
    carboxamido = Chem.MolFromSmarts("N-C(=O)")
    if not mol.HasSubstructMatch(carboxamido):
        return False, "Missing a carboxamido group (N-C(=O) fragment)"

    return True, "Molecule has a penam core with two methyl substituents at position 2, a carboxylate at position 3, and a carboxamido group at position 6."

# Example usage (if run as a script)
if __name__ == "__main__":
    # Test with one of the provided penicillin SMILES (e.g., Penicillin K)
    test_smiles = "CCCCCCCC(=O)N[C@H]1[C@H]2SC(C)(C)[C@@H](N2C1=O)C(O)=O"
    is_pen, reason = is_penicillin(test_smiles)
    print("Is penicillin?", is_pen)
    print("Reason:", reason)