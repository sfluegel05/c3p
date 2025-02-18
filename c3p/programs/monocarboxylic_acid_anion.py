"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
#!/usr/bin/env python3
"""
Classifies: Monocarboxylic acid anion
Definition: A carboxylic acid anion formed when the carboxyl group of a monocarboxylic acid is deprotonated.
A “true” monocarboxylic acid anion should have exactly one [C(=O)[O-]] group, an overall formal charge of –1,
and no other atoms carry nonzero formal charge. Furthermore, the carboxyl carbon is expected to be trigonal,
bonded to exactly two oxygens (one in a C=O carbonyl with formal charge 0 and one terminal O (i.e. deprotonated)
with formal charge –1) and one R‐group. Note that while these rules capture typical examples, many borderline
cases exist and expert review may be necessary.
"""

from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    The algorithm works in several steps:
      • Parse the SMILES.
      • Search for exactly one carboxylate group matching the SMARTS "[#6X3](=O)[O-]".
      • Check that the overall formal charge of the molecule is exactly –1.
      • Verify that only one atom carries a negative charge.
      • Confirm that the carboxyl carbon has exactly three neighbors, two of which are oxygens.
      • Determine (from the bond orders) which oxygen is double‐bonded (carbonyl; charge 0)
        and which is single‐bonded (deprotonated; charge –1) and additionally require that the 
        deprotonated oxygen is terminal (degree 1).
      • (Extra check) Confirm that the remaining non‐oxygen neighbor (the “R‐group”) is unique.
    If any step fails, a (False, reason) pair is returned.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple; if True the molecule is classified as a monocarboxylic acid anion,
                     otherwise False with an explanation.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for a deprotonated carboxylate group:
    # It should be a trigonal carbon with one double bonded oxygen and one oxygen with a -1 charge.
    carboxylate_smarts = "[#6X3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    if carboxylate_pattern is None:
        return False, "Error generating SMARTS for carboxylate group"
    
    # Look for carboxylate substructure matches.
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(matches) == 0:
        return False, "No deprotonated carboxyl (carboxylate) group found"
    if len(matches) > 1:
        return False, f"Found {len(matches)} carboxylate groups; molecule is not a monocarboxylic acid anion"
    
    # Check overall formal charge.
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != -1:
        return False, f"Expected overall charge of -1 for a monocarboxylate anion, found charge = {total_charge}"
    
    # Ensure that exactly one atom carries a negative formal charge.
    neg_charged_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0]
    if len(neg_charged_atoms) != 1:
        return False, ("More than one atom carries a negative charge; expected a single deprotonated oxygen "
                       "for a monocarboxylic acid anion")
    
    # The carboxylate match returns three atom indices: (carboxyl carbon, oxygen, oxygen)
    match = matches[0]
    if len(match) != 3:
        return False, "Unexpected match size for carboxylate group; expected 3 atoms"
    carboxyl_idx, oxy1_idx, oxy2_idx = match
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    oxy1_atom = mol.GetAtomWithIdx(oxy1_idx)
    oxy2_atom = mol.GetAtomWithIdx(oxy2_idx)
    
    # Extra check: the carboxyl carbon should have exactly 3 neighbors.
    if carboxyl_atom.GetDegree() != 3:
        return False, f"Carboxyl carbon (atom idx {carboxyl_idx}) does not have 3 neighbors; connectivity is unexpected"
    
    # Among its neighbors, exactly 2 must be oxygens.
    o_neighbors = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(o_neighbors) != 2:
        return False, "Carboxyl carbon does not have exactly 2 oxygen neighbors"
    
    # Identify which oxygen is bonded via a double bond (carbonyl, expected formal charge 0)
    # versus a single bond (the deprotonated oxygen, expected -1 charge).
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
    
    # Check that the deprotonated oxygen is terminal (degree should be 1).
    if deprot_oxygen.GetDegree() != 1:
        return False, "Deprotonated oxygen is not terminal as expected"
    
    # Check formal charges on the oxygens.
    if carbonyl_oxygen.GetFormalCharge() != 0:
        return False, "Carbonyl oxygen does not have the expected formal charge of 0"
    if deprot_oxygen.GetFormalCharge() != -1:
        return False, "Deprotonated oxygen does not carry a formal charge of -1"
    
    # Confirm that the only negatively charged atom in the molecule is the deprotonated oxygen.
    if neg_charged_atoms[0].GetIdx() != deprot_oxygen.GetIdx():
        return False, ("An unexpected atom carries the negative charge (should be solely on the deprotonated oxygen)")
    
    # Extra check: the carboxyl carbon should have exactly one non-oxygen neighbor (the R‐group).
    non_o_neighbors = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() != 8]
    if len(non_o_neighbors) != 1:
        return False, "Carboxyl carbon does not have exactly one non-oxygen neighbor (R‐group expected)"
    # (Optional: one might further inspect the R-group to see if it fits the expected profile for a simple monocarboxylic acid.)
    
    # Finally, verify that no other atoms (outside the matched carboxylate group) possess nonzero formal charge.
    carboxylate_atom_ids = set(match)
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 and atom.GetIdx() not in carboxylate_atom_ids:
            return False, ("Additional charged atoms found outside the carboxylate group, "
                           "which is inconsistent with a simple monocarboxylic acid anion")
    
    # If we reached here, the molecule passed all our tests.
    return True, ("Contains a single deprotonated carboxylate group (with expected connectivity and charge distribution), "
                  "consistent with a monocarboxylic acid anion. (Note: borderline cases may require expert review.)")

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Example: xanthine-8-carboxylate should be a true positive.
    test_smiles = "[O-]C(=O)c1nc2[nH]c(=O)[nH]c(=O)c2[nH]1"
    result, reason = is_monocarboxylic_acid_anion(test_smiles)
    print(result, reason)