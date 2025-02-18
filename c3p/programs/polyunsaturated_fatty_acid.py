"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: Polyunsaturated Fatty Acid
Definition: A fatty acid with a terminal –COOH group that is attached to a single, 
            linear, non‐aromatic aliphatic chain and containing more than one 
            non-aromatic carbon–carbon double bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    The molecule must contain a terminal carboxylic acid group (the acid carbon is
    attached to exactly one carbon) and a long, linear, non-aromatic fatty acyl chain 
    with more than one non-aromatic C=C double bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule qualifies as a polyunsaturated fatty acid, False otherwise.
        str: Explanation for classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a carboxylic acid group. We use SMARTS that matches a –COOH unit.
    # [CX3](=O)[O;H] matches a carbon with one double bond O and an OH.
    ca_smarts = "[CX3](=O)[O;H]"
    ca_group = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    if not ca_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # Check for a terminal acid: the acid carbon should be attached to exactly one carbon.
    terminal_acid_found = False
    acid_atom = None  # carbonyl carbon of the –COOH
    acyl_chain_start = None  # carbon attached to the acid that starts the chain
    for match in ca_matches:
        # In our pattern, the first atom (index 0 in the match) is the acid (carbonyl) carbon.
        acid_c = mol.GetAtomWithIdx(match[0])
        # Count the carbon neighbors of the acid carbon (ignore oxygens).
        carbon_neighbors = [nbr for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_acid_found = True
            acid_atom = acid_c
            acyl_chain_start = carbon_neighbors[0]
            break
    if not terminal_acid_found or acyl_chain_start is None:
        return False, "Carboxylic acid group not terminal; not a typical fatty acid"
    
    # Now we want to “grow” the fatty acyl chain starting from acyl_chain_start.
    # We require that the chain consists of non-aromatic carbons that are not in rings.
    def longest_chain_path(atom, coming_from, visited):
        """
        Recursively finds the longest linear chain (as a list of atom indices) starting from 'atom'
        while avoiding going back to the 'coming_from' atom. Only non-aromatic, non-ring carbon 
        atoms are allowed.
        """
        best_path = [atom.GetIdx()]
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == coming_from.GetIdx():
                continue  # do not go back to previous atom
            # Only consider carbon atoms that are non-aromatic and not in any ring.
            if nbr.GetAtomicNum() != 6 or nbr.GetIsAromatic() or nbr.IsInRing():
                continue
            if nbr.GetIdx() in visited:
                continue
            new_visited = visited | {nbr.GetIdx()}
            candidate_path = longest_chain_path(nbr, atom, new_visited)
            candidate_path = [atom.GetIdx()] + candidate_path
            if len(candidate_path) > len(best_path):
                best_path = candidate_path
        return best_path

    # Start the search from acyl_chain_start; do not let it go back into the acid.
    visited = {acid_atom.GetIdx(), acyl_chain_start.GetIdx()}
    chain_path = longest_chain_path(acyl_chain_start, acid_atom, visited)
    chain_length = len(chain_path)
    
    # Set a threshold for what is considered a "long chain".
    MIN_CHAIN_LENGTH = 8
    if chain_length < MIN_CHAIN_LENGTH:
        return False, f"Fatty acid chain length only {chain_length} carbons; too short to be a typical fatty acid"
    
    # Count the number of non-aromatic C=C double bonds along the acyl chain.
    # We examine bonds between consecutive atoms in the longest chain path.
    double_bond_count = 0
    for i in range(len(chain_path) - 1):
        bond = mol.GetBondBetweenAtoms(chain_path[i], chain_path[i+1])
        if bond:
            if bond.GetBondType() == Chem.BondType.DOUBLE and (not bond.GetIsAromatic()):
                double_bond_count += 1

    if double_bond_count <= 1:
        return False, f"Found {double_bond_count} non-aromatic C=C bond(s) in chain; need more than one to qualify as polyunsaturated"
    
    return True, f"Contains a terminal carboxylic acid group attached to a chain of {chain_length} carbons and {double_bond_count} non-aromatic double bonds; qualifies as a polyunsaturated fatty acid"

# Example usage:
if __name__ == "__main__":
    # Try one of the provided examples: resolvin D6.
    test_smiles = "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O"
    result, reason = is_polyunsaturated_fatty_acid(test_smiles)
    print(result, reason)