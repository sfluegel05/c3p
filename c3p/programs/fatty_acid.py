"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: fatty acid
Definition: “Any acyclic aliphatic monocarboxylic acid. Natural fatty acids commonly have a carbon chain of 4 to 28 carbons 
(usually unbranched and even‐numbered). By extension, the term is sometimes used to embrace all acyclic aliphatic carboxylic acids.”
This improved classifier:
  - Requires that a free (undissociated) carboxylic acid group –C(=O)[OH] is present.
  - Ensures that the acid carbon is “clean” (besides its two oxygens it should only have one carbon neighbour).
  - From that unique neighbour (the “α‐carbon”) it DFS–searches over all simple paths (only through carbon atoms not in rings 
    and that do not show extra carbonyl bonds) and takes the longest path – thus allowing for branched chains.
  - The total count (acid carbon + chain length) must be between 4 and 28.
  
If a candidate is found, returns True plus an explanation of the chain length.
Otherwise returns False with a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if the given SMILES represents a fatty acid as defined above.
    Rules:
      - Valid molecule.
      - No nitrogen atoms (to avoid peptides/amino acid derivatives).
      - Contains a free (undissociated) carboxyl group –C(=O)[OH].
      - The carboxyl carbon (the [CX3] in [CX3](=O)[OX2H]) must have exactly one carbon neighbour,
        so that the acid group attaches to a single candidate aliphatic chain.
      - From that carbon neighbour, the algorithm recursively finds the longest acyclic carbon-only path 
        (ignoring atoms in rings and paths where an extra carbonyl functionality is present).
      - The total chain “length” is defined as 1 (the acid carbon) plus the number of carbons along the longest path.
        Fatty acids must have a total chain length between 4 and 28 (inclusive).
      
    Args:
      smiles (str): SMILES string for the molecule.
      
    Returns:
      bool: True if the molecule is classified as a fatty acid, False otherwise.
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Reject molecules with nitrogen atoms (signals peptides/amino acid derivatives)
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Molecule contains nitrogen atoms indicating peptide or amino acid derivatives"
    
    # SMARTS to find free carboxylic acid group –C(=O)[OH]
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No free carboxyl (–C(=O)[OH]) group detected"
    
    # Helper: for a candidate atom (except the acid carbon) check if it displays an extra carbonyl functionality.
    # Extra C=O bonds (double bonds to oxygen) outside the acid group mark extra functionalities.
    def has_extra_carbonyl(atom, exclude_bonds=set()):
        for bond in atom.GetBonds():
            if bond.GetIdx() in exclude_bonds:
                continue
            # Check if the bond is double and the other atom is oxygen.
            if bond.GetBondTypeAsDouble() == 2.0:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    return True
        return False

    # DFS function: from a starting carbon atom, explore all simple paths on carbon atoms.
    # Only follow atoms not in rings and that do not show extra C=O bonds.
    # Returns the maximum number of carbons in a valid path (counting the starting atom as 1).
    def dfs_chain(atom, visited):
        # If this atom shows an extra carbonyl (apart from acid functionality) then disqualify this branch.
        if has_extra_carbonyl(atom):
            return 0
        max_length = 1  # count current atom
        for nbr in atom.GetNeighbors():
            # Only follow carbon neighbours
            if nbr.GetAtomicNum() != 6:
                continue
            # Do not follow atoms in rings (the candidate chain must be acyclic)
            if nbr.IsInRing():
                continue
            if nbr.GetIdx() in visited:
                continue
            # Explore further along this branch.
            new_visited = visited.copy()
            new_visited.add(nbr.GetIdx())
            branch_length = 1 + dfs_chain(nbr, new_visited)
            if branch_length > max_length:
                max_length = branch_length
        return max_length

    # For each acid candidate, try to find a valid candidate chain.
    for match in acid_matches:
        acid_carbon_idx = match[0]  # the [CX3] atom in the acid pattern
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
        
        # Identify all carbon neighbours of the acid carbon.
        # In a proper free fatty acid the carboxyl carbon should only be attached to its two oxygens and a single carbon.
        carbon_neigh = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neigh) != 1:
            continue  # this acid candidate has extra substituents on the acid carbon
        
        alpha = carbon_neigh[0]
        if alpha.IsInRing():
            continue  # reject if the chain starts in a ring
        
        # Start DFS from the alpha carbon.
        visited = {acid_carbon_idx, alpha.GetIdx()}
        chain_candidate_length = dfs_chain(alpha, visited)
        # Total chain length = acid carbon (1) + chain candidate length.
        total_chain_length = 1 + chain_candidate_length
        
        if 4 <= total_chain_length <= 28:
            reason = (f"Found free acid group (–C(=O)[OH]) with a candidate contiguous aliphatic chain of "
                      f"{total_chain_length} carbons (including the acid carbon).")
            return True, reason

    return False, "No candidate free acid group is attached to an acceptable contiguous aliphatic chain (4–28 carbons)"

# Example test cases (executing when running this script)
if __name__ == '__main__':
    test_set = [
        # True positives
        ("OC(=O)CCCCCCCCCCC=C", "12-tridecenoic acid"),
        ("C(C(O)=O)C/C=C\\C[C@H](\\C=C\\C=C/C=C/C=C/[C@H]([C@@H](C/C=C\\CC)O)O)O", "aspirin-triggered resolvin D2"),
        ("O[C@H](CCC)C/C=C\\C/C=C\\CCCCCCCC(O)=O", "Avenoleic acid"),
        ("OC(=O)C=CC=CC=CCCCCCCCCCCCCCCCCC", "Tetracosatrienoic acid"),
        ("C(CCCCCCC[C@H]([C@H](CCCCCCCC)O)O)(=O)O", "(9R,10S)-dihydroxyoctadecanoic acid"),
        ("CCCCCC\\C=C/CCCCCC(O)=O", "cis-tetradec-7-enoic acid"),
        ("OC(=O)CCCCCCCCCCC/C=C/CC", "13-hexadecenoic acid"),
        ("CCCCCC[C@@H](O)\\C=C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O", "15(R)-HETE"),
        ("[H]\\C(C)=C\\C(O)=O", "isocrotonic acid"),
        # False positives (expected to be rejected)
        ("CCCCCC/C=C\\C/C=C\\CCCCCCCCCCCCC(O)=O", "Prenateic acid"),
        ("OC(CCCCCCCCCC/C=C\\CCC)=O", "Criegeenic acid"),
        ("C[C@@H]1O[C@@H](OCCCCCCCCCCCCCCCCCCCC[C@@H](O)CC(O)=O)[C@H](O)C[C@H]1O", "bhos#42"),
        ("C[C@@H]1O[C@@H](OCCCCCC(O)=O)[C@H](O)C[C@H]1O", "oscr#12"),
        # False negatives (should be accepted as fatty acids)
        ("OC1C(C(C(=O)C1)C/C=C\\C/C=C\\CC)/C=C/C(O)C/C=C\\CCC(O)=O", "7-hydroxy-D4-neuroprostane"),
        ("CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O", "juvenile hormone I acid"),
        ("O=C(O)[C@H]([C@@H](C(=O)O)C)CCCCCCCC", "Sphaeric acid"),
        ("CC(C)C[C@@H](O)C(O)=O", "(R)-2-hydroxy-4-methylpentanoic acid"),
    ]
    
    for smi, name in test_set:
        result, explanation = is_fatty_acid(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nExplanation: {explanation}\n")