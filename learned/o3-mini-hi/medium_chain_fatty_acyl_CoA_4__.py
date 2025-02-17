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
  5. The acyl chain attached at the thioester – that is, from the carbonyl carbon towards the R‐group (excluding the sizeable CoA core) – must be “medium‐chain” meaning that, when counting only linear (acyclic) carbon atoms (including the carbonyl carbon) and not including any parts of the CoA core, the number is between 6 and 12.
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
    thioester_smarts = "[C;!R](=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) bond found"
    # Take the first match (a tuple of atom indices); the first atom is the carbonyl carbon and the next is the sulfur.
    carbonyl_idx, sulfur_idx = thioester_matches[0]
    
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
    # We gather atom indices belonging to the adenine and pantetheine fragments.
    excluded_atoms = set()
    for match in mol.GetSubstructMatches(adenine_pattern):
        excluded_atoms.update(match)
    for match in mol.GetSubstructMatches(pantetheine_pattern):
        excluded_atoms.update(match)
    # (Note: we intentionally do not remove these atoms so that we keep connectivity; rather we ignore them during chain search.)
    
    # --- (5) Extract the acyl chain ---
    # In a fatty acyl-CoA, the acyl group is R-CO–, where the CO is the carbonyl we matched.
    # Its neighboring carbon (other than the sulfur, which is part of the CoA) belongs to the acyl chain.
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    acyl_neighbors = []
    for nbr in carbonyl_atom.GetNeighbors():
        # We want a carbon neighbor that is not (likely) part of the CoA core.
        if nbr.GetSymbol() == "C" and nbr.GetIdx() not in excluded_atoms:
            acyl_neighbors.append(nbr)
    if not acyl_neighbors:
        return False, "No acyl chain attached to the thioester bond"
    
    # Define a DFS that follows only linear (acyclic) carbon paths.
    def dfs_chain_length(atom, visited):
        # Count the current atom (each is a carbon)
        max_len = 1  
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() != "C":
                continue
            if nbr.IsInRing():
                continue
            if nbr.GetIdx() in visited:
                continue
            # Do not go into excluded CoA atoms.
            if nbr.GetIdx() in excluded_atoms:
                continue
            # Continue depth-first.
            new_visited = visited.copy()
            new_visited.add(nbr.GetIdx())
            branch_len = 1 + dfs_chain_length(nbr, new_visited)
            if branch_len > max_len:
                max_len = branch_len
        return max_len

    # For each allowed neighbor from the carbonyl, compute the chain length.
    candidate_chain_length = 0
    for acyl_atom in acyl_neighbors:
        # Only follow acyclic paths (ignore if the neighbor is in a ring even though unlikely for fatty acids)
        if acyl_atom.IsInRing():
            continue
        # The chain length includes the carbonyl carbon (hence add 1); then follow the neighbor chain.
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
    # Running some example SMILES strings (from the provided outcomes)
    test_smiles_list = [
        # True positives:
        "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)C/C=C\\CCCCCC)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O",  # cis-dec-3-enoyl-CoA(4-)
        "CCCCCCCCCC(C)=OSCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12",  # undecanoyl-CoA(4-)
        # A couple false positives/negatives can be tested here as well.
    ]
    for smi in test_smiles_list:
        result, reason = is_medium_chain_fatty_acyl_CoA_4__(smi)
        print(result, reason)