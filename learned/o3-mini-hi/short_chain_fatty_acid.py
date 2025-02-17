"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: Short‐chain fatty acid.
Definition:
  "An aliphatic monocarboxylic acid with a chain length of less than C6.
   If any non‐hydrocarbon substituent is present, the compound is not normally regarded as a short‐chain fatty acid."
Our implementation requires that:
  – the molecule is acyclic
  – exactly one carboxyl group is present (using the SMARTS "[CX3](=O)[OX2H]")
  – the COOH carbon is terminal (it is attached to exactly one carbon atom)
  – all oxygens are found in that carboxyl group (i.e. the total oxygen count equals 2)
  – the total number of carbon atoms (i.e. the “chain length” including the acid carbon) is at least 3 and at most 6.
  
These additional cuts help to reject both acids that are “too small” (e.g. formic acid, acetic acid)
and those that, although they pass simpler structural tests, are not normally considered short‐chain fatty acids.
"""

from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as a short-chain fatty acid.
    
    A qualifying short-chain fatty acid is defined as:
      (a) an aliphatic (acyclic) molecule,
      (b) containing exactly one carboxylic acid group ([CX3](=O)[OX2H]),
      (c) in which the carboxyl carbon is terminal (attached to exactly one other carbon),
      (d) with no extra heteroatoms outside the carboxyl group (i.e. the only oxygens are those in COOH),
      (e) and the total number of carbons is between 3 and 6 (inclusive).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a short-chain fatty acid, False otherwise.
        str: A message describing the reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Must be fully acyclic
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Molecule contains rings; expected an acyclic (aliphatic) acid."
    
    # Identify the carboxyl group using SMARTS.
    ca_group_smarts = "[CX3](=O)[OX2H]"
    ca_group = Chem.MolFromSmarts(ca_group_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    if not ca_matches:
        return False, "No carboxylic acid functional group found."
    if len(ca_matches) != 1:
        return False, f"Expected exactly one carboxyl group, found {len(ca_matches)}."
    
    # By our SMARTS the first atom of the match is the carboxyl carbon.
    carboxyl_idx = ca_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    
    # Check that the carboxyl carbon is terminal in the carbon-only graph:
    # Count how many neighboring atoms (by bonds) that are carbon atoms.
    carbon_neighbors = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, ("The carboxyl carbon is not terminal; "
                       f"expected exactly one carbon neighbor, found {len(carbon_neighbors)}.")
    
    # Count the total number of oxygen atoms in the molecule.
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:
        return False, f"Expected exactly 2 oxygens (from the carboxyl group), found {oxygen_count}."
    
    # Count the total number of carbon atoms in the molecule.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Also impose a lower bound so that trivial acids like formic (1 C) or acetic acid (2 C) are excluded.
    if carbon_count < 3:
        return False, f"Total carbon count is {carbon_count}, which is too low for a fatty acid."
    if carbon_count > 6:
        return False, f"Total carbon count is {carbon_count}, which exceeds allowed short-chain length (max 6)."
    
    # Now check that no extra (non-hydrocarbon) substituents occur outside the carboxyl group.
    # (We allow only C and H in the remainder; note that hydrogen atoms are usually implicit.)
    ca_atoms = set(ca_matches[0])
    for atom in mol.GetAtoms():
        if atom.GetIdx() in ca_atoms:
            continue  # atoms in the carboxyl group are allowed (carboxyl oxygens)
        # Only allow carbon atoms.
        if atom.GetAtomicNum() != 6:
            return False, (f"Found a non-hydrocarbon substituent: atom {atom.GetSymbol()} "
                           "outside the carboxyl group is not allowed.")
    
    return True, ("Molecule contains one carboxyl group with a terminal acid carbon, is acyclic, "
                  "has 3–6 carbons total, and no extraneous substituents.")

# Example usage:
if __name__ == '__main__':
    # A list of SMILES strings taken from the test outcomes.
    test_cases = [
        ("CC(C)(C)C(O)=O", "pivalic acid"),  # correct: 5 carbons
        ("CCC(O)=O", "propionic acid"),       # correct: 3 carbons
        ("C=CCC(C)C(=O)O", "2-Methyl-4-pentenoic acid"), # should be false positive
        ("OC(=O)C(CC)=CC", "2-ethyl-2-butenoic acid"),   # false positive expected
        ("OC(=O)/C=C(\\CC)/C", "3-methyl-2Z-pentenoic acid"),  # false positive expected
        ("[H]C(O)=O", "formic acid"),         # false: too few carbons
        ("CC(O)=O", "acetic acid"),           # false: too few carbons
        ("CCCC(O)=O", "butyric acid"),         # correct: 4 carbons
        ("C[CH]C(C)C(O)=O", "2-methylbutyric acid"),  # correct if 5 carbons (example)
        ("C/C=C/C(O)=O", "cis-pent-2-enoic acid"),     # correct: 5 carbons even with unsaturation
    ]
    
    for s, name in test_cases:
        valid, reason = is_short_chain_fatty_acid(s)
        print(f"SMILES: {s} NAME: {name} -> {valid}; Reason: {reason}")