"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: medium-chain fatty acyl-CoA(4-)
Definition: An acyl-CoA oxoanion that results from deprotonation of the phosphate and diphosphate groups 
of any medium-chain fatty acyl-CoA; major species at pH 7.3.
Criteria:
  1. SMILES must be valid.
  2. The molecule must contain a thioester bond (a non‐ring carbonyl bound to sulfur).
  3. It must contain two characteristic substructures of Coenzyme A:
       - an adenine (nucleobase) ring,
       - and a pantetheine-like fragment (SCCNC(=O)CCNC(=O)).
  4. The overall formal charge must be –4.
  5. The acyl chain attached at the thioester – that is, from the carbonyl carbon towards the R‐group 
     (excluding the sizeable CoA core) – must be “medium‐chain” meaning that, when counting only linear 
     (acyclic) carbon atoms (including the carbonyl carbon) and not including any parts of the CoA core, 
     the number is between 6 and 12.
If any criterion fails, the function returns False with a diagnostic explanation.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule matches the medium-chain fatty acyl-CoA(4-) pattern, else False.
        str: Explanation for the decision.
    """
    # Parse SMILES; if invalid, abort.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- (1) Check for the thioester bond ---
    # We look for a non‐ring carbonyl (C not in a ring) doubly bonded to O and bonded to a sulfur.
    # The SMARTS "[C;!R](=O)S" captures three atoms: carbon, oxygen, and sulfur.
    thioester_smarts = "[C;!R](=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) bond found"
    # thioester_matches[0] returns a tuple of three atom indices. We take the first (carbonyl carbon)
    # and the last (sulfur) for our purpose.
    match = thioester_matches[0]
    if len(match) < 3:
        return False, "Unexpected match size for thioester bond"
    carbonyl_idx = match[0]
    sulfur_idx = match[-1]
    
    # --- (2) Check for CoA core fragments ---
    # Adenine (nucleobase) pattern.
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine (nucleobase) fragment found"
    # Pantetheine fragment.
    pantetheine_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine fragment (SCCNC(=O)CCNC(=O)) found"
    
    # --- (3) Check overall formal charge ---
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -4:
        return False, f"Net charge is {net_charge}, but expected -4 for this class"
    
    # --- (4) Exclude CoA core atoms for the acyl chain extraction ---
    excluded_atoms = set()
    for match in mol.GetSubstructMatches(adenine_pattern):
        excluded_atoms.update(match)
    for match in mol.GetSubstructMatches(pantetheine_pattern):
        excluded_atoms.update(match)
    # We intentionally do not remove these atoms from the molecule; we simply ignore them during the chain search.
    
    # --- (5) Extract the acyl chain ---
    # The fatty acyl chain is the fragment attached to the thioester bond.
    # From the thioester, starting at the matched carbonyl carbon, get its neighbors excluding the sulfur.
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    acyl_neighbors = []
    for nbr in carbonyl_atom.GetNeighbors():
        # Exclude the sulfur that is part of the thioester bond.
        if nbr.GetIdx() == sulfur_idx:
            continue
        if nbr.GetSymbol() == "C" and nbr.GetIdx() not in excluded_atoms:
            acyl_neighbors.append(nbr)
    if not acyl_neighbors:
        return False, "No acyl chain attached to the thioester bond"
    
    # We'll perform a depth-first search (DFS) to find the maximum length of a linear (acyclic) carbon chain.
    def dfs_chain_length(atom, visited):
        max_len = 1  # count the current carbon atom in the chain
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() != "C":
                continue
            if nbr.IsInRing():
                continue  # ignore atoms in rings
            if nbr.GetIdx() in visited:
                continue
            if nbr.GetIdx() in excluded_atoms:
                continue  # skip atoms that belong to the CoA core
            new_visited = visited.copy()
            new_visited.add(nbr.GetIdx())
            branch_len = 1 + dfs_chain_length(nbr, new_visited)
            if branch_len > max_len:
                max_len = branch_len
        return max_len

    candidate_chain_length = 0
    # For each allowed neighbor from the carbonyl, compute chain length.
    for acyl_atom in acyl_neighbors:
        if acyl_atom.IsInRing():
            continue
        chain_length = 1 + dfs_chain_length(acyl_atom, visited={carbonyl_idx, acyl_atom.GetIdx()})
        if chain_length > candidate_chain_length:
            candidate_chain_length = chain_length

    if candidate_chain_length == 0:
        return False, "Could not extract an acyl chain from the thioester bond"
    
    # --- (6) Check that the acyl chain length is in the medium-chain range (6-12 carbons) ---
    if candidate_chain_length < 6 or candidate_chain_length > 12:
        return False, (f"Acyl chain length is {candidate_chain_length} carbons, "
                       "not in medium-chain range (6-12)")
    
    return True, (f"Matches acyl-CoA pattern with medium chain length of {candidate_chain_length} carbons "
                  "and net charge -4")

# Example usage:
if __name__ == "__main__":
    # We include a couple of example SMILES strings from the provided outcomes.
    test_smiles_list = [
        # Example: cis-dec-3-enoyl-CoA(4-)
        "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)C/C=C\\CCCCCC)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O",
        # Example: undecanoyl-CoA(4-) (this example may not exactly match the SMARTS for CoA core fragments,
        # but it is provided as a test case).
        "CCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    ]
    for smi in test_smiles_list:
        result, reason = is_medium_chain_fatty_acyl_CoA_4__(smi)
        print(result, reason)