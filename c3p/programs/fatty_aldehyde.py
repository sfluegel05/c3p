"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: Fatty Aldehyde
Definition: An aldehyde arising from reduction of a fatty acid’s –COOH group,
with a long acyclic chain that is mainly carbon/hydrogen and a terminal aldehyde group.
Improvements over the previous version:
  • Only allows elements C, H, and O.
  • Requires an acyclic structure.
  • Demands at least 6 carbon atoms.
  • Requires the heavy atoms (all atoms except H) be predominantly carbon (fraction >= 0.75).
  • Requires exactly one overall carbonyl group ([CX3]=[OX1]) in the molecule.
  • Within that carbonyl, the aldehyde carbon must be terminal (attached to exactly one carbon).
"""

from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule qualifies as a fatty aldehyde.

    Steps:
      1. Parse the SMILES string.
      2. Confirm that only the allowed elements (C, H, O) are present.
      3. Check that the molecule is acyclic.
      4. Verify that there are at least 6 carbon atoms.
      5. Confirm that the heavy atom (non-H) carbon fraction is at least 0.75.
      6. Verify that exactly one carbonyl group ([CX3]=[OX1]) is present.
      7. From the carbonyl matches, ensure that a terminal aldehyde exists.
         Terminal means that the aldehyde carbon is not in a ring and is attached (ignoring the =O)
         to exactly one other carbon.
    Args:
      smiles (str): SMILES string of the molecule.
    Returns:
      (bool, str): A tuple where the boolean indicates acceptance and the string explains the decision.
    """

    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- Check allowed elements: only C (6), H (1), and O (8) ---
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element: {atom.GetSymbol()}"

    # --- Check that the molecule is acyclic ---
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s) which is not typical for a fatty aldehyde"

    # --- Count the number of carbon atoms ---
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(carbon_atoms)
    if n_carbons < 6:
        return False, f"Not enough carbon atoms ({n_carbons} found; need at least 6) for a fatty chain"

    # --- Check heavy atom composition (non-H atoms) ---
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not heavy_atoms:
        return False, "No heavy atoms found"
    carbon_fraction = n_carbons / len(heavy_atoms)
    if carbon_fraction < 0.75:
        return False, f"Molecule has a low carbon fraction ({carbon_fraction:.2f}); likely contains additional functionalities"

    # --- Check that there is exactly one carbonyl group overall ---
    # SMARTS for a carbonyl group (C=O) of any type.
    carbonyl_smarts = "[CX3]=[OX1]"
    carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) != 1:
        return False, f"Expected exactly one carbonyl group; found {len(carbonyl_matches)}"

    # --- Check for a terminal aldehyde group ---
    # Define SMARTS for an aldehyde: a carbon with exactly one H and double-bonded to O.
    aldehyde_smarts = "[CX3H1](=O)"
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    terminal_aldehyde_count = 0
    for match in aldehyde_matches:
        # The first index is the aldehyde carbon.
        aldehyde_c = mol.GetAtomWithIdx(match[0])
        if aldehyde_c.IsInRing():
            continue  # Exclude if the carbonyl carbon is in a ring.
        # Count only carbon neighbors (ignore the double-bonded oxygen).
        carbon_neighbors = [nbr for nbr in aldehyde_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_aldehyde_count += 1

    if terminal_aldehyde_count != 1:
        return False, f"Expected exactly one terminal aldehyde group; found {terminal_aldehyde_count}"

    return True, "Molecule qualifies as a fatty aldehyde: acyclic, mainly C/H, with a long chain and a single terminal aldehyde group"

# Uncomment the following lines to run a few example test cases:
# test_smiles = [
#     "[H]C(=CC=O)C(O)CCCCC",    # 4-hydroxynon-2-enal (expected True)
#     "O=CCCCCCCCCC/C=C\\CCCCCCCC",  # 11Z-Eicosenal (expected True)
#     "CCCCCCCCCCCC=O",         # dodecanal (expected True)
#     "O=C/C(=C\\CCCCCC/C=C\\C/C=C\\CCCCC)/CCCCCCCCCCCCC",  # (2Z,10Z,13Z)-2-tridecylnonadeca-2,10,13-trienal (expected True)
# ]
# for s in test_smiles:
#     res, reason = is_fatty_aldehyde(s)
#     print(f"SMILES: {s}\nResult: {res}, Reason: {reason}\n")