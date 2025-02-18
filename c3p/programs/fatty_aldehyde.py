"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: Fatty Aldehyde
Definition: An aldehyde arising from the reduction of the carboxylic acid group of its 
corresponding fatty acid, having a carbonyl group at one end of a long (typically acyclic) carbon chain.
Improvements:
  - Only allows molecules made of C, H, and O.
  - Accepts acyclic molecules.
  - Requires the molecule to have at least six carbons.
  - Demands a high fraction of carbons among heavy atoms.
  - Requires exactly one terminal aldehyde group (a [CX3H1](=O) where the carbon is attached to exactly one carbon).
"""

from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    
    Checks performed:
      1. Valid SMILES parsing.
      2. The molecule is composed only of allowed elements (C, H, O).
      3. The molecule is acyclic (contains no rings).
      4. The molecule contains at least 6 carbon atoms.
      5. The heavy atoms should be dominated by carbons (carbon fraction > 0.8).
      6. The molecule contains exactly one terminal aldehyde group.
         Here, a terminal aldehyde is defined as a [CX3H1](=O) group where the aldehyde carbon:
           - Is not in a ring.
           - Is connected (excluding the double-bonded oxygen) to exactly one other carbon atom.
      
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if molecule is classified as a fatty aldehyde, False otherwise.
        str: Explanation for the decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check that no disallowed elements (only allow C, H, O) ---
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element: {atom.GetSymbol()}"
    
    # --- Check acyclicity ---
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s) which is not typical for a fatty aldehyde"
    
    # --- Count total number of carbon atoms ---
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(carbon_atoms)
    if n_carbons < 6:
        return False, f"Not enough carbon atoms ({n_carbons} found; need at least 6) for a fatty chain"
    
    # --- Check heavy atom composition: count atoms with atomic number > 1 ---
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    n_heavy = len(heavy_atoms)
    if n_heavy == 0:
        return False, "No heavy atoms found"
    carbon_fraction = n_carbons / n_heavy
    if carbon_fraction < 0.8:
        return False, f"Molecule has a low carbon fraction ({carbon_fraction:.2f}); likely contains additional functionalities"
    
    # --- Check for terminal aldehyde group ---
    # Define SMARTS for aldehyde: a carbon with one hydrogen and double-bonded to oxygen.
    aldehyde_smarts = "[CX3H1](=O)"
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not matches:
        return False, "No aldehyde group found"
    
    # We require exactly one terminal aldehyde.
    terminal_aldehyde_count = 0
    for match in matches:
        # match[0] is the aldehyde carbon.
        aldehyde_c = mol.GetAtomWithIdx(match[0])
        if aldehyde_c.IsInRing():
            continue  # aldehyde carbon must not be in a ring
        
        # Count neighbors that are carbons (excluding the oxygen double-bonded)
        carbon_neighbors = [nbr for nbr in aldehyde_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue  # terminal aldehyde should be attached to exactly one carbon
        
        # Passed this check: count it as a valid terminal aldehyde 
        terminal_aldehyde_count += 1

    if terminal_aldehyde_count != 1:
        return False, f"Expected exactly one terminal aldehyde group; found {terminal_aldehyde_count}"
    
    # Passed all the checks
    return True, "Molecule qualifies as a fatty aldehyde: acyclic, mainly carbon/hydrogen, with a long chain and a terminal aldehyde group"

# Uncomment to run some example test‐cases:
# test_smiles = [
#     "[H]C(=CC=O)C(O)CCCCC",   # 4-hydroxynon-2-enal, expected true
#     "O=CCCCCCCCCC/C=C\\CCCCCCCC",  # 11Z-Eicosenal, expected true
#     "CCCCCCCCCCCC=O",         # dodecanal, expected true
#     "[H]C(=O)C([H])=C([H])CCCC",  # hept-2-enal, now accepted (7 carbons)
#     "O=CC=CCCC=CCC",         # 2,6-nonadienal, expected true
# ]
# for s in test_smiles:
#     res, reason = is_fatty_aldehyde(s)
#     print(f"SMILES: {s}\nResult: {res}, Reason: {reason}\n")