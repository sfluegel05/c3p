"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
#!/usr/bin/env python3
"""
Classifies: Monocarboxylic acid anion
Definition: A carboxylic acid anion formed when the carboxyl group of a monocarboxylic acid is deprotonated.
That is, the molecule should have exactly one deprotonated carboxyl (CO2–) group with the proper connectivity,
an overall formal charge of -1, and no other atoms in the molecule should be charged.
Extra connectivity checks are performed to ensure that the carboxylate carbon is attached to exactly two oxygen atoms
(with one oxygen double‐bonded and neutral, and one single‐bonded bearing the -1 charge) and one R‐group.
"""

from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    
    Criteria:
      1. The SMILES must parse.
      2. It must contain exactly one deprotonated carboxylate group. We locate that group via the SMARTS:
           "[#6X3](=O)[O-]"
         which requires that the carboxyl carbon (#6) is trigonal (X3), has a double bond to one oxygen,
         and a single bond to an oxygen carrying a -1 charge.
      3. The overall formal charge of the molecule must be -1.
      4. No atom outside the matched carboxylate group may carry a nonzero formal charge.
      5. Extra connectivity is checked: the carboxyl carbon should have exactly 3 neighbors,
         exactly two of which are oxygens. Furthermore, one of these oxygens should be connected via a double bond
         (carbonyl oxygen, formal charge 0) and the other via a single bond (deprotonated oxygen, formal charge -1).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): Tuple where the boolean indicates if the molecule is classified as a monocarboxylic acid anion,
                     and the string provides an explanation.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a deprotonated carboxylate group.
    carboxylate_smarts = "[#6X3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    if carboxylate_pattern is None:
        return False, "Error creating SMARTS for carboxylate group"

    # Find all matches for the carboxylate pattern.
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(matches) == 0:
        return False, "No deprotonated carboxyl (carboxylate) group found"
    if len(matches) > 1:
        return False, f"Found {len(matches)} carboxylate groups; molecule is not a monocarboxylic acid anion"

    # Check overall formal charge.
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != -1:
        return False, f"Expected overall charge of -1 for a monocarboxylate anion, found charge = {total_charge}"

    # Identify the atoms involved in the carboxylate group.
    # We expect the SMARTS match to return a tuple: (carboxyl carbon, carbonyl oxygen, deprotonated oxygen)
    match = matches[0]
    if len(match) != 3:
        return False, "Unexpected match size for carboxylate group; expected 3 atoms"
    carboxyl_idx, oxy1_idx, oxy2_idx = match
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    oxy1_atom = mol.GetAtomWithIdx(oxy1_idx)
    oxy2_atom = mol.GetAtomWithIdx(oxy2_idx)

    # Extra connectivity check: carboxyl carbon should have exactly 3 neighbors.
    if carboxyl_atom.GetDegree() != 3:
        return False, (f"Carboxyl carbon (atom idx {carboxyl_idx}) does not have 3 neighbors; "
                       "unexpected connectivity for a simple acid group")
    
    # Among its neighbors, exactly two must be oxygens.
    o_neighbors = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(o_neighbors) != 2:
        return False, "The carboxyl carbon does not have exactly 2 oxygen neighbors"
    
    # Check bond types to the oxygen neighbors:
    # One oxygen should be connected via a double bond (the carbonyl oxygen, expected to have formal charge 0)
    # and the other via a single bond (the deprotonated oxygen, expected to have formal charge -1).
    bond1 = mol.GetBondBetweenAtoms(carboxyl_idx, oxy1_idx)
    bond2 = mol.GetBondBetweenAtoms(carboxyl_idx, oxy2_idx)
    if bond1 is None or bond2 is None:
        return False, "Bond information missing for carboxylate group"

    # Identify the double vs single bond.
    # RDKit represents bond types as Chem.BondType.DOUBLE and .SINGLE.
    if bond1.GetBondType() == Chem.BondType.DOUBLE and bond2.GetBondType() == Chem.BondType.SINGLE:
        carbonyl_oxygen = oxy1_atom
        deprot_oxygen = oxy2_atom
    elif bond2.GetBondType() == Chem.BondType.DOUBLE and bond1.GetBondType() == Chem.BondType.SINGLE:
        carbonyl_oxygen = oxy2_atom
        deprot_oxygen = oxy1_atom
    else:
        return False, "The bonds from carboxyl carbon to oxygen atoms do not display the expected single/double pattern"

    # Check the formal charges on the oxygens: carbonyl oxygen should be neutral (0) and the deprotonated oxygen -1.
    if carbonyl_oxygen.GetFormalCharge() != 0:
        return False, "Carbonyl oxygen has an unexpected formal charge (expected 0)"
    if deprot_oxygen.GetFormalCharge() != -1:
        return False, "Deprotonated oxygen does not have a formal charge of -1"

    # Now check that no other atoms (outside the carboxylate group) carry a nonzero formal charge.
    carboxylate_atom_ids = set(match)
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 and atom.GetIdx() not in carboxylate_atom_ids:
            return False, ("Additional charged atoms found outside the deprotonated carboxyl group, "
                           "indicating presence of other ionizable groups")
    
    # Passed all tests.
    return True, ("The molecule contains exactly one deprotonated carboxyl group (with proper connectivity), "
                  "an overall charge of -1, and no other charged atoms, consistent with a monocarboxylic acid anion")


# Example usage (can be removed if used as a module)
if __name__ == "__main__":
    # Example: xanthine-8-carboxylate should be classified as a monocarboxylic acid anion.
    test_smiles = "[O-]C(=O)c1nc2[nH]c(=O)[nH]c(=O)c2[nH]1"
    result, reason = is_monocarboxylic_acid_anion(test_smiles)
    print(result, reason)