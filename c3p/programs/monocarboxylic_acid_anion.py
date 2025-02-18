"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
#!/usr/bin/env python3
"""
Classifies: Monocarboxylic acid anion
Definition: A carboxylic acid anion formed when the carboxyl group of a monocarboxylic acid is deprotonated.
That is, the molecule should have exactly one deprotonated carboxyl group (a carboxylate group)
and an overall formal charge of -1. In addition, no other atoms should be charged and
the carboxylate group should have the typical connectivity of a simple acid (a carbonyl carbon attached to exactly two oxygens and one R-group).
"""

from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on a SMILES string.
    
    The criteria are:
      1. The molecule must parse.
      2. It must contain exactly one deprotonated carboxyl group. We use a SMARTS pattern that forces
         the carboxyl carbon to be a carbon atom (atomic number 6) with trigonal geometry:
           "[#6X3](=O)[O-]"
      3. The overall formal charge of the molecule must be -1.
      4. No atom outside the matched carboxylate group should have a nonzero formal charge.
      5. The carboxylate carbon should have exactly three neighbors (two oxygens from the carboxylate and one R‐group).
         This extra connectivity check helps avoid classifying more complex acyl substructures.
    
    Args:
        smiles (str): The SMILES string representing the molecule.

    Returns:
         (bool, str): True plus a message if the molecule is classified as a monocarboxylic acid anion,
                      otherwise False and an explanatory message.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for deprotonated carboxylate; require the acidic carbon to be a carbon (#6)
    carboxylate_smarts = "[#6X3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    if carboxylate_pattern is None:
        return False, "Error creating SMARTS pattern for carboxylate group"
    
    # Find all matches of the carboxylate pattern
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(matches) == 0:
        return False, "No deprotonated carboxyl (carboxylate) group found"
    if len(matches) > 1:
        return False, f"Found {len(matches)} carboxylate groups; molecule is not a monocarboxylic acid anion"
    
    # Check overall formal charge: must be exactly -1
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != -1:
        return False, f"Expected overall charge of -1 for a monocarboxylate anion, found charge = {total_charge}"

    # Get the indices from the single carboxylate match.
    # We assume the SMARTS matches in order:
    #   index 0: carboxyl carbon, index 1: carbonyl oxygen, index 2: deprotonated oxygen.
    carboxylate_atoms = set(matches[0])
    carboxyl_c = matches[0][0]
    # Extra connectivity check for the carboxyl carbon:
    # It should have exactly three neighbors: the two oxygens in the carboxylate and one additional R-group atom.
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_c)
    if carboxyl_atom.GetDegree() != 3:
        return False, f"Carboxyl carbon (atom idx {carboxyl_c}) does not have 3 neighbors, unexpected connectivity for a simple acid group"
    
    # (Optional) Check that none of the non-carboxylate atoms have nonzero formal charge.
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 and atom.GetIdx() not in carboxylate_atoms:
            return False, ("Additional charged atoms found outside the deprotonated carboxyl group, "
                           "indicating the presence of other ionizable groups")
    
    # Passed all tests; we classify this as a monocarboxylic acid anion.
    return True, ("The molecule contains exactly one deprotonated carboxyl group (with proper connectivity) "
                  "and an overall charge of -1 with no other charged atoms, consistent with a monocarboxylic acid anion")

# Example usage (can be removed if used as a module):
if __name__ == "__main__":
    test_smiles = "[O-]C(=O)c1nc2[nH]c(=O)[nH]c(=O)c2[nH]1"  # xanthine-8-carboxylate example
    result, reason = is_monocarboxylic_acid_anion(test_smiles)
    print(result, reason)