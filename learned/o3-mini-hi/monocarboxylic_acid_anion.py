"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
#!/usr/bin/env python3
"""
Classifies: Monocarboxylic acid anion
Definition: A carboxylic acid anion formed when the carboxyl group of a monocarboxylic acid is deprotonated.
The molecule must have exactly one deprotonated carboxylate group (with proper connectivity), an overall
formal charge of -1, and no other atom should carry a nonzero formal charge.
Extra connectivity checks enforce that the carboxyl carbon is trigonal (three neighbors), 
with exactly two oxygens (one double-bonded carbonyl oxygen having 0 charge and one single-bonded oxygen with -1 charge)
and exactly one “R‐group” (typically carbon).
"""

from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.

    This revised algorithm checks not only that there is exactly one carboxylate group matching the
    SMARTS "[#6X3](=O)[O-]", and that the overall formal charge equals -1, but also that:
      • the deprotonated oxygen is the only atom bearing a negative charge
      • the carboxyl carbon is connected to exactly two oxygens (with the proper single/double bond pattern)
        and one non-oxygen (the R‐group).
    
    Caveat: Some structures may have resonance or extra substituents that cause ambiguity and thus
    extra filtering (or manual inspection) may be necessary.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): Tuple where the boolean indicates whether the molecule is classified as a 
                     monocarboxylic acid anion, and the string explains the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Create a SMARTS pattern for a deprotonated carboxylate group: a trigonal carboxyl carbon
    # double bonded to one oxygen and single-bonded to an oxygen carrying a -1 charge.
    carboxylate_smarts = "[#6X3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    if carboxylate_pattern is None:
        return False, "Error generating SMARTS for carboxylate group"
    
    # Find all matches to the carboxylate pattern.
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(matches) == 0:
        return False, "No deprotonated carboxyl (carboxylate) group found"
    if len(matches) > 1:
        return False, f"Found {len(matches)} carboxylate groups; molecule is not a monocarboxylic acid anion"
    
    # Check that the overall formal charge of the molecule is exactly -1.
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != -1:
        return False, f"Expected overall charge of -1 for a monocarboxylate anion, found charge = {total_charge}"
    
    # Enforce that there is exactly one atom in the molecule with a negative formal charge.
    neg_charged_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0]
    if len(neg_charged_atoms) != 1:
        return False, ("More than one atom carries a negative charge; "
                       "expected a single deprotonated oxygen for a monocarboxylic acid anion")
    
    # The SMARTS match returns three atom indices: (carboxyl carbon, oxygen, oxygen)
    match = matches[0]
    if len(match) != 3:
        return False, "Unexpected match size for carboxylate group; expected 3 atoms"
    carboxyl_idx, oxy1_idx, oxy2_idx = match
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    oxy1_atom = mol.GetAtomWithIdx(oxy1_idx)
    oxy2_atom = mol.GetAtomWithIdx(oxy2_idx)
    
    # Extra connectivity check: the carboxyl carbon should have exactly 3 neighbors.
    if carboxyl_atom.GetDegree() != 3:
        return False, f"Carboxyl carbon (atom idx {carboxyl_idx}) does not have 3 neighbors; connectivity is unexpected"
    
    # Among its neighbors, exactly 2 must be oxygens.
    o_neighbors = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(o_neighbors) != 2:
        return False, "Carboxyl carbon does not have exactly 2 oxygen neighbors"
    
    # Identify which oxygen is bonded via a double bond (the carbonyl oxygen, expected charge 0)
    # versus a single bond (the deprotonated oxygen with -1 charge).
    bond1 = mol.GetBondBetweenAtoms(carboxyl_idx, oxy1_idx)
    bond2 = mol.GetBondBetweenAtoms(carboxyl_idx, oxy2_idx)
    if bond1 is None or bond2 is None:
        return False, "Missing bond information in carboxylate group"
    
    if bond1.GetBondType() == Chem.BondType.DOUBLE and bond2.GetBondType() == Chem.BondType.SINGLE:
        carbonyl_oxygen = oxy1_atom
        deprot_oxygen = oxy2_atom
    elif bond2.GetBondType() == Chem.BondType.DOUBLE and bond1.GetBondType() == Chem.BondType.SINGLE:
        carbonyl_oxygen = oxy2_atom
        deprot_oxygen = oxy1_atom
    else:
        return False, "Bonds from carboxyl carbon to oxygens do not match the expected double/single pattern"
    
    # Check formal charges on the oxygen atoms.
    if carbonyl_oxygen.GetFormalCharge() != 0:
        return False, "Carbonyl oxygen does not have the expected formal charge of 0"
    if deprot_oxygen.GetFormalCharge() != -1:
        return False, "Deprotonated oxygen does not carry a formal charge of -1"
    
    # Confirm that the only negatively charged atom in the molecule is the deprotonated oxygen.
    if neg_charged_atoms[0].GetIdx() != deprot_oxygen.GetIdx():
        return False, "An unexpected atom carries the negative charge (should be solely on the deprotonated oxygen)"
    
    # Finally, verify that no other atoms (outside the carboxylate group) possess nonzero formal charge.
    carboxylate_atom_ids = set(match)
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 and atom.GetIdx() not in carboxylate_atom_ids:
            return False, ("Additional charged atoms found outside the carboxylate group, "
                           "which is inconsistent with a simple monocarboxylic acid anion")
    
    # All tests passed.
    return True, ("The molecule contains a single deprotonated carboxylate group with proper connectivity, "
                  "an overall charge of -1, and no other charged atoms, as expected for a monocarboxylic acid anion.")

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Test with an example (xanthine-8-carboxylate, which should be classified as True)
    test_smiles = "[O-]C(=O)c1nc2[nH]c(=O)[nH]c(=O)c2[nH]1"
    result, reason = is_monocarboxylic_acid_anion(test_smiles)
    print(result, reason)