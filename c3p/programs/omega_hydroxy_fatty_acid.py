"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: Omega-hydroxy fatty acid
Definition: A naturally-occurring straight-chain fatty acid that is composed only of C, H, and O,
with a terminal carboxyl group (–COOH) at position 1 and a hydroxyl group (–OH) at the omega (last) carbon.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is an omega-hydroxy fatty acid.
    The criteria required are:
      1. Only the elements C, H, and O are present.
      2. The molecule is acyclic.
      3. There is a carboxylic acid (–COOH) group, identified using SMARTS "[CX3](=O)[OX2H]".
      4. The carbon atom of the –COOH is one terminal atom of an unbranched (linear) carbon chain.
      5. The chain is identified by enumerating all simple paths in the carbon subgraph from the
         carboxyl carbon; the longest such path is taken as the fatty acid’s main chain.
      6. The omega (last) carbon (the other terminus of that chain) should bear a hydroxyl (-OH) group.
      
    Args:
        smiles (str): SMILES string for the molecule.
  
    Returns:
        bool: True if the molecule meets the criteria for an omega-hydroxy fatty acid, False otherwise.
        str: Explanation for the determination.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure only elements C, H, and O are present.
    allowed = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed:
            return False, "Molecule contains atoms other than C, H, and O"

    # Exclude cyclic molecules (straight-chain requirement)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic; expected a straight-chain structure"

    # Add explicit hydrogens to help recognize -OH groups.
    mol = Chem.AddHs(mol)

    # Find the carboxylic acid group using SMARTS.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid (-COOH) group found in the molecule"

    # Use the first match as the carboxyl group. We assume the first atom (index 0 in the match tuple)
    # is the carboxyl carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)

    # In a fatty acid, the acid carbon (the one in -COOH) should be attached to only one carbon.
    acid_carbons = [nbr.GetIdx() for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(acid_carbons) != 1:
        return False, "Carboxylic acid group is not at a terminal position (acid carbon is bonded to >1 carbon)"

    # Build the subgraph of all carbon atoms.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Create a dictionary mapping carbon index -> list of neighboring carbon indices.
    carbon_graph = {}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        carbon_graph[idx] = neighbors

    # To capture the main fatty acid chain, we enumerate all simple paths (acyclic)
    # starting at the acid carbon. Each simple path is a candidate chain.
    # In a straight-chain fatty acid the acid carbon will be at one terminus.
    all_paths = []
    def dfs(current, path, visited):
        # Mark the current node as visited.
        visited.add(current)
        path.append(current)
        # If current is a terminal in the carbon_graph (only one neighbor) and it is not the starting acid carbon,
        # record this path.
        if (current != acid_carbon_idx) and (len(carbon_graph[current]) == 1):
            all_paths.append(path.copy())
        else:
            # Continue DFS on neighbors that have not been visited.
            for nbr in carbon_graph[current]:
                if nbr not in visited:
                    dfs(nbr, path, visited)
        path.pop()
        visited.remove(current)

    dfs(acid_carbon_idx, [], set())
    
    if not all_paths:
        return False, "No continuous carbon chain found from the carboxyl group"

    # Choose the longest path as the main chain.
    longest = max(all_paths, key=lambda p: len(p))
    
    # Ensure that the acid carbon is actually at one end of this longest path.
    if acid_carbon_idx != longest[0] and acid_carbon_idx != longest[-1]:
        return False, "Carboxylic acid group is not located at a terminal position of the carbon chain"
    
    # Identify the omega (last) carbon as the terminus opposite the acid carbon.
    omega_candidate_idx = longest[-1] if acid_carbon_idx == longest[0] else longest[0]
    omega_candidate = mol.GetAtomWithIdx(omega_candidate_idx)

    # Check that the omega candidate carries an -OH substituent.
    has_omega_oh = False
    for nbr in omega_candidate.GetNeighbors():
        if nbr.GetAtomicNum() == 8:  # oxygen atom
            # Check that the oxygen is part of a hydroxyl (has at least one hydrogen attached)
            if any(sub.GetAtomicNum() == 1 for sub in nbr.GetNeighbors()):
                has_omega_oh = True
                break

    if not has_omega_oh:
        return False, "Terminal (omega) carbon lacks a hydroxyl (-OH) substituent"

    return True, "Molecule is a straight-chain fatty acid with a COOH group at C1 and an -OH at the omega position"

    
# Example test (remove or comment out these test calls as needed):
if __name__ == '__main__':
    # Example: (2E)-18-hydroxyoctadec-2-enoic acid SMILES.
    test_smiles = "OCCCCCCCCCCCCCCC\\C=C\\C(O)=O"
    result, reason = is_omega_hydroxy_fatty_acid(test_smiles)
    print(result, reason)