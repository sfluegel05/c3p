"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid – any straight-chain C18 polyunsaturated fatty acid
with exactly 2 non-aromatic C=C double bonds and a terminal carboxylic acid group.
"""

from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string)
    is an octadecadienoic acid. The enforced definition is:
      • The molecule possesses a carboxylic acid group at one terminus.
      • Its main (longest) carbon chain -- forced to start at the acid group --
        is a straight (unbranched) chain of exactly 18 carbon atoms.
      • Along that chain there are exactly 2 non-aromatic C=C bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple. True with a positive message if it is an octadecadienoic acid,
                     otherwise False with an explanatory reason.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find a terminal carboxylic acid group.
    # The SMARTS below matches a carboxyl carbon: [CX3](=O)[OX2H1]
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Missing terminal carboxylic acid functionality"
    
    # Build an undirected graph connecting only carbon atoms.
    carbon_graph = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # count only carbons
            idx = atom.GetIdx()
            # Only add other carbons as neighbors.
            carbon_graph[idx] = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    
    # To ensure the chain is “straight” we check that every carbon in the candidate chain
    # has no extra carbon neighbors other than those in the linear sequence.
    # Terminal carbons in a straight chain must have exactly one carbon neighbor,
    # and internal carbons exactly two.
    def is_straight_chain(chain):
        for i, idx in enumerate(chain):
            # Get all carbon neighbors (from the full molecule’s carbon graph).
            c_neighbors = set(carbon_graph.get(idx, []))
            # Neighbors that belong to the chain.
            in_chain_neighbors = set()
            if i > 0:
                in_chain_neighbors.add(chain[i-1])
            if i < len(chain)-1:
                in_chain_neighbors.add(chain[i+1])
            # In a straight chain no additional carbon neighbor should exist.
            if c_neighbors != in_chain_neighbors:
                return False
        return True

    # For a fatty acid the acid carbon should be one terminus.
    # We collect all the acid carbon indices from the SMARTS match.
    acid_carbon_set = {match[0] for match in acid_matches if match}
    
    # We now search for a straight-chain path of exactly 18 carbons that begins with one
    # of these acid carbons. (The choice of starting at the acid group reduces ambiguity.)
    found_chain = None  # store a candidate chain if found

    def dfs(path):
        nonlocal found_chain
        if found_chain is not None:  # solution already found
            return
        if len(path) == 18:
            # Candidate chain found; now check if it is unbranched (i.e. "straight")
            if is_straight_chain(path):
                found_chain = path[:]  # record a copy of the path
            return
        last = path[-1]
        for nbr in carbon_graph.get(last, []):
            # Do not revisit nodes already in the chain
            if nbr in path:
                continue
            dfs(path + [nbr])
            if found_chain is not None:
                return

    # Run DFS starting from each acid carbon.
    for acid_c in acid_carbon_set:
        dfs([acid_c])
        if found_chain is not None:
            break

    if found_chain is None:
        return False, "No straight chain of 18 carbons found starting from the carboxylic acid group"
    
    # Now count the non-aromatic C=C bonds along the candidate chain.
    double_bond_count = 0
    for i in range(len(found_chain) - 1):
        bond = mol.GetBondBetweenAtoms(found_chain[i], found_chain[i+1])
        if bond is None:
            return False, f"Missing bond between atoms {found_chain[i]} and {found_chain[i+1]}"
        # Count only double bonds that are non-aromatic.
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and not bond.GetIsAromatic():
            double_bond_count += 1

    if double_bond_count != 2:
        return False, f"Found {double_bond_count} non-aromatic C=C bonds along the straight chain; exactly 2 required"
    
    return True, ("Molecule is a straight-chain C18 fatty acid with exactly 2 non-aromatic C=C bonds "
                  "and a terminal carboxylic acid group")

# When run as a script, a few examples can be tested.
if __name__ == "__main__":
    test_cases = [
        # True positives (should be classified as octadecadienoic acid)
        ("CCCCCC\\C=C/C=C/[C@H](O)CCCCCCCC(O)=O", "9(R)-HODE"),
        ("OC(=O)CCCCCCC\\C=C/C=C\\CCCCCC", "9Z,11Z-octadecadienoic acid"),
        ("CC/C=C\\C/C=C\\CC(C(CCCCCCCC(O)=O)O)O", "9,10-DiHODE"),
        ("OC(CCCCCCCC(O)=O)C(O)/C=C/C(O)C/C=C\\CC", "(11E,15Z)-9,10,13-trihydroxyoctadeca-11,15-dienoic acid"),
        ("O[C@H](CCC)C/C=C\\C/C=C\\CCCCCCCC(O)=O", "Avenoleic acid"),
        ("O[C@H](C/C=C\\CCCCC)\\C=C\\CCCCCCC(O)=O", "10R-HODE"),
        # False positives and negatives would be similarly tested.
        ("C(=C\\C/C=C\\CCCCCO)CCCCCCCC(=O)O", "18-hydroxylinoleic acid (should be rejected)"),
    ]

    for sm, name in test_cases:
        res, reason = is_octadecadienoic_acid(sm)
        status = "CORRECT" if res else "REJECTED"
        print(f"SMILES: {sm}\nNAME: {name}\nResult: {status}\nReason: {reason}\n")