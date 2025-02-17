"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: Polyunsaturated Fatty Acid
Definition: Any fatty acid (as judged by having a terminal carboxylic acid group that 
            is attached to a long, mostly aliphatic chain) containing more than one 
            non-aromatic carbon-carbon double bond.
This improved version attempts to limit false positives by ensuring that the acid group 
is at the terminus of a long chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    The molecule must contain a terminal carboxylic acid group (i.e. the acid carbon is
    attached to a single, long, non-aromatic chain) and more than one non-aromatic
    carbon-carbon double bond.
    
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
    
    # First, check for a carboxylic acid group.
    # We use a SMARTS that matches a terminal –COOH group:
    # This pattern ensures that the acid carbon (CX3) is double-bonded to an oxygen and bonded to an –OH.
    # Later, we will check if it is attached to exactly one non-oxygen carbon.
    ca_smarts = "[CX3](=O)[O;H]"
    ca_group = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    if not ca_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # Look for a terminal carboxylic acid: the acid carbon should have exactly one neighboring carbon
    terminal_acid_found = False
    acyl_chain_start = None  # the carbon that attaches to the carboxyl group
    for match in ca_matches:
        # In our pattern, the first atom is the carbonyl carbon.
        acid_c = mol.GetAtomWithIdx(match[0])
        # Count how many carbon neighbors (ignoring oxygen atoms) acid_c has.
        carbon_neighbors = [nbr for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_acid_found = True
            acyl_chain_start = carbon_neighbors[0]
            break
    if not terminal_acid_found or acyl_chain_start is None:
        return False, "Carboxylic acid group not terminal; not a typical fatty acid"

    # Define a helper function to compute the longest contiguous carbon chain starting
    # from the given atom by traversing along non-aromatic carbon atoms. We avoid revisiting atoms.
    def get_longest_chain(atom, coming_from, visited):
        longest = 1  # count the current atom
        for nbr in atom.GetNeighbors():
            # avoid going backwards
            if nbr.GetIdx() == coming_from.GetIdx():
                continue
            # Only consider carbons that are not aromatic.
            if nbr.GetAtomicNum() != 6 or nbr.GetIsAromatic():
                continue
            if nbr.GetIdx() in visited:
                continue
            new_visited = visited | {nbr.GetIdx()}
            branch_len = 1 + get_longest_chain(nbr, atom, new_visited)
            if branch_len > longest:
                longest = branch_len
        return longest

    # Obtain the longest chain length starting from the atom directly attached to the terminal acid.
    chain_length = get_longest_chain(acyl_chain_start, coming_from=mol.GetAtomWithIdx(acyl_chain_start.GetNeighbors()[0].GetIdx()), visited={acyl_chain_start.GetIdx()})
    # Set a threshold for what is considered a "long chain" for a fatty acid (typical fatty acids are generally long).
    MIN_CHAIN_LENGTH = 8
    if chain_length < MIN_CHAIN_LENGTH:
        return False, f"Fatty acid chain length only {chain_length} carbons; too short to be a typical fatty acid"

    # Next, count the number of non-aromatic carbon-carbon double bonds.
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            if bond.GetIsAromatic():
                continue
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Count only if both atoms are carbons.
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                double_bond_count += 1

    if double_bond_count <= 1:
        return False, f"Found {double_bond_count} non-aromatic carbon-carbon double bond(s); need more than one to qualify as polyunsaturated"

    return True, f"Contains a terminal carboxylic acid group attached to a chain of {chain_length} carbons and {double_bond_count} non-aromatic double bonds; qualifies as a polyunsaturated fatty acid"

# Example usage:
if __name__ == "__main__":
    test_smiles = "C(C(O)=O)C/C=C\\C[C@H](\\C=C\\C=C/C=C/C=C/[C@H]([C@@H](C/C=C\\CC)O)O)O"  # resolvin D6
    result, reason = is_polyunsaturated_fatty_acid(test_smiles)
    print(result, reason)