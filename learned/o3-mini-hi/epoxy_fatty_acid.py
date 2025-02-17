"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: Epoxy fatty acid
Definition: A heterocyclic fatty acid containing an epoxide ring as part of its structure.
Improvement: Instead of rejecting extra rings outright (which led to false negatives),
we now require that the molecule contains a terminal carboxylic acid group,
a long acyclic (mostly chain-like) region, and exactly one 3-membered epoxide ring.
Additional rings are allowed provided the acyclic chain is sufficiently long 
and the overall ratio of acyclic carbons is high.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    Improved criteria:
      - Contains a terminal carboxylic acid group (C(=O)[O;H]) with the acid carbon attached to only one other carbon.
      - Has a sufficiently long aliphatic chain. Here we trace, from the neighbor of the acid carbon,
        the longest path along acyclic carbon atoms (ignoring bonds that are in rings) and require at least 12 carbons.
      - Contains exactly one epoxide ring (a 3-membered ring with one oxygen and two carbons, matched by [CX2]1[OX2][CX2]1).
      - The overall carbon skeleton should be predominantly acyclic. We compute the ratio of acyclic carbon atoms vs.
        all carbon atoms and require it be at least 0.5.
      - Molecule’s molecular weight should be above 200 Da (typical for fatty acids).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as an epoxy fatty acid, False otherwise.
        str: A reason message explaining the classification.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find a terminal carboxylic acid group.
    # SMARTS for carboxylic acid: carbon double-bonded to oxygen and single-bonded to an OH.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found to signify a fatty acid"
    
    terminal_acid_found = False
    acid_carbon_idx = None
    acid_neighbor_idx = None
    for match in acid_matches:
        # match[0] is the acid carbon.
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # Look for neighboring carbons (ignore O neighbors).
        carbon_neighbors = [nbr.GetIdx() for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_acid_found = True
            acid_carbon_idx = acid_carbon.GetIdx()
            acid_neighbor_idx = carbon_neighbors[0]
            break
    if not terminal_acid_found:
        return False, "Carboxylic acid group is found but not terminal"
    
    # Check basic fatty acid properties: number of carbon atoms and molecular weight.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 12:
        return False, f"Too few carbon atoms ({total_carbons}) to be a fatty acid"
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.2f} Da) too low for a fatty acid"
    
    # Compute the ratio of acyclic carbons.
    acyclic_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing())
    acyclic_ratio = acyclic_carbons / total_carbons if total_carbons > 0 else 0
    if acyclic_ratio < 0.5:
        return False, f"Low acyclic carbon ratio ({acyclic_ratio:.2f}); structure too cyclic to be a fatty acid"
    
    # Helper: from a starting atom, traverse along bonds between carbons that are not in rings.
    visited = set()
    def dfs(atom_idx):
        visited.add(atom_idx)
        max_length = 1  # count the starting atom
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            # Only follow if neighbor is carbon, not in a ring, and not already visited.
            if nbr.GetAtomicNum() == 6 and not nbr.IsInRing() and nbr.GetIdx() not in visited:
                length = 1 + dfs(nbr.GetIdx())
                if length > max_length:
                    max_length = length
        visited.remove(atom_idx)
        return max_length
    
    # Begin DFS from the neighbor of the terminal acid carbon.
    chain_length = dfs(acid_neighbor_idx)
    if chain_length < 12:
        return False, f"Longest acyclic carbon chain from the acid group is too short (length: {chain_length})"
    
    # Identify epoxide rings: 3-membered rings with exactly 1 oxygen and 2 carbons.
    # We use SMARTS that matches an epoxide: two sp3 carbons and one oxygen in a ring.
    epoxide_pattern = Chem.MolFromSmarts("[CX2]1[OX2][CX2]1")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    # To avoid double‐counting, we only count unique sets (sorted tuples).
    epoxide_sets = set(tuple(sorted(match)) for match in epoxide_matches)
    if len(epoxide_sets) != 1:
        return False, f"Expected exactly one epoxide ring, found {len(epoxide_sets)}"
    
    return True, "Contains a terminal carboxylic acid group, a long acyclic chain, and one epoxide ring indicative of an epoxy fatty acid"

# Example usage:
# Uncomment the following lines to test:
# test_smiles = "C(CCC/C=C\\C[C@@H]1/C(/O1)=C/C=C\\C/C=C\\CCCCC)(=O)O"  # (5Z,8R,9Z,11Z,14Z)-8,9-epoxyicosatetraenoic acid
# result, reason = is_epoxy_fatty_acid(test_smiles)
# print(result, reason)