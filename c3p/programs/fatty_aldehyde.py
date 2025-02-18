"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: Fatty Aldehyde
Definition: An aldehyde formally arising from reduction of the carboxylic acid group 
of its corresponding fatty acid, having a carbonyl group at one end of the carbon chain.
We expect a fatty aldehyde to be acyclic, with most heavy atoms being carbons (and O’s from a few 
hydroxyls), and a single terminal aldehyde group.
"""

from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    
    Checks performed:
      1. The SMILES must be valid.
      2. The molecule must contain at least 8 carbons.
      3. The molecule must be acyclic (no rings), as fatty aldehydes are typically linear.
      4. The molecule must contain at least one aldehyde group [CX3H1](=O) that is terminal:
         - The aldehyde carbon should not be in a ring.
         - It must have exactly one carbon neighbor (attached to the aliphatic chain).
      5. The overall heavy-atom composition should be dominated by C (and allowed O) –
         if the fraction of carbons among heavy atoms is low, the molecule likely contains other 
         functionalities and is mis‐classified.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if molecule is classified as a fatty aldehyde, False otherwise.
        str: Explanation for the decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total number of carbon atoms
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(carbons)
    if n_carbons < 8:
        return False, f"Not enough carbon atoms ({n_carbons} found; need at least 8) for a fatty chain"
    
    # Check if the molecule is acyclic (no rings); fatty aldehydes are expected to be linear
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s) which is not typical for a fatty aldehyde"
    
    # Compute heavy-atom count (exclude hydrogens) and carbon fraction.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    n_heavy = len(heavy_atoms)
    if n_heavy == 0:
        return False, "No heavy atoms found"
    carbon_fraction = n_carbons / n_heavy
    # In typical fatty aldehydes most heavy atoms are carbons (and one or two oxygens)
    if carbon_fraction < 0.75:
        return False, "Molecule has a low carbon fraction (likely contains additional functionalities)"
    
    # Define aldehyde SMARTS: a carbon with one hydrogen and double-bonded to oxygen.
    aldehyde_smarts = "[CX3H1](=O)"
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"
    
    # Check each aldehyde match to ensure it is terminal (i.e., attached to exactly one carbon)
    terminal_aldehyde_found = False
    for match in aldehyde_matches:
        # match[0] is the aldehyde carbon (the one bearing the hydrogen).
        aldehyde_c = mol.GetAtomWithIdx(match[0])
        if aldehyde_c.IsInRing():
            continue  # must not be in a ring
        
        # Find non-oxygen neighbors (the C=O oxygen is not used to define the chain).
        neighbor_carbons = []
        for nb in aldehyde_c.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                continue
            if nb.GetAtomicNum() == 6:
                neighbor_carbons.append(nb)
        # Terminal aldehyde should have exactly one carbon neighbor.
        if len(neighbor_carbons) != 1:
            continue
        
        # Additionally, check that the neighbor carbon (start of the chain) is not in a ring.
        chain_carbon = neighbor_carbons[0]
        if chain_carbon.IsInRing():
            continue
        
        # Passed the terminal aldehyde checks: we consider this aldehyde as valid.
        terminal_aldehyde_found = True
        break
    
    if not terminal_aldehyde_found:
        return False, "No valid terminal aldehyde group found; aldehyde may be internal or part of a ring"
    
    return True, "Molecule qualifies as a fatty aldehyde: has a long acyclic carbon chain with a terminal aldehyde group"

# Example test‐cases (uncomment to run):
# test_smiles = [
#     "[H]C(=CC=O)C(O)CCCCC",  # 4-hydroxynon-2-enal, true positive
#     "O=CCCCCCCCCC/C=C\\CCCCCCCC",  # 11Z-Eicosenal, true positive
#     "CCCCCCCCCCCC=O",  # dodecanal, true positive
#     "[H]C(=O)C([H])=C([H])CCCC",  # hept-2-enal; false negative if chain has only 7 C's
# ]
# for s in test_smiles:
#     res, reason = is_fatty_aldehyde(s)
#     print(f"SMILES: {s} -> {res}: {reason}")