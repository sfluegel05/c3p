"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: medium-chain fatty acyl-CoA(4-)
Definition: An acyl-CoA oxoanion that results from deprotonation of the phosphate and diphosphate groups 
of any medium-chain fatty acyl-CoA; major species at pH 7.3.
The criteria here are:
  1. SMILES must be valid.
  2. The molecule must contain a thioester bond (a non‐ring carbonyl bound to sulfur).
  3. It must contain two characteristic substructures of Coenzyme A:
       - an adenine (nucleobase) ring,
       - and a pantetheine-like fragment (SCCNC(=O)CCNC(=O)) present in most acyl‐CoAs.
  4. The overall formal charge must be -4.
  5. The acyl chain attached at the thioester must be “medium‐chain” in that,
     when counting only linear (acyclic) carbon atoms (starting at the carbonyl carbon), 
     the chain length (including the carbonyl carbon) is between 6 and 12.
     
If any of these criteria fail, the function returns False with a diagnostic explanation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is medium-chain fatty acyl-CoA(4-), else False.
        str: Explanation for the decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the thioester bond.
    # Use a SMARTS that demands the carbonyl is not in a ring: this helps avoid false positives.
    # [C;!R](=O)S matches a nonring carbonyl carbon doubly bonded to oxygen and linked to a sulfur.
    thioester_smarts = "[C;!R](=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) bond found"
    
    # 2. Check for CoA fragments.
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")  # adenine-type ring(s)
    pantetheine_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine (nucleobase) fragment found"
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine fragment (SCCNC(=O)CCNC(=O)) found"
    
    # 3. Check overall formal charge.
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -4:
        return False, f"Net charge is {net_charge}, but expected -4 for this class"

    # 4. Remove the CoA core fragments to better isolate the acyl chain.
    # We remove the adenine and pantetheine fragments, hoping that the acyl chain still remains attached
    # to the thioester carbonyl.
    mol_clean = Chem.DeleteSubstructs(mol, adenine_pattern)
    mol_clean = Chem.DeleteSubstructs(mol_clean, pantetheine_pattern)
    # Sometimes deleting substructures leaves dummy atoms (e.g. [*]) that we remove.
    mol_clean = Chem.RemoveHs(Chem.RemoveHs(mol_clean))
    
    # 5. Extract the acyl chain.
    # Choose one thioester match. In the match tuple, the first atom is the carbonyl carbon.
    carbonyl_idx = thioester_matches[0][0]
    # Try to locate this carbonyl in the cleaned molecule.
    # (We compare based on atom mapping; if the deletion failed to remove the carbonyl, we use original.)
    try:
        carbonyl_atom_clean = mol_clean.GetAtomWithIdx(carbonyl_idx)
    except Exception:
        carbonyl_atom_clean = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Among the neighbors of the carbonyl, choose the one that is a carbon (and not the sulfur)
    acyl_neighbors = [nbr for nbr in carbonyl_atom_clean.GetNeighbors() if nbr.GetSymbol() == "C"]
    if not acyl_neighbors:
        return False, "No acyl chain attached to the thioester bond"
    
    # Use a DFS routine that only walks over carbon atoms not in rings.
    # (This is our heuristic for the “linear” chain.)
    def dfs_chain_length(atom, visited):
        cur_idx = atom.GetIdx()
        visited.add(cur_idx)
        max_len = 1  # count current carbon
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() != "C":
                continue
            if nbr.IsInRing():
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            path_len = 1 + dfs_chain_length(nbr, visited.copy())
            if path_len > max_len:
                max_len = path_len
        return max_len

    candidate_chain_length = 0
    # We treat each acyl neighbor as the start of a potential chain.
    for acyl_atom in acyl_neighbors:
        if acyl_atom.IsInRing():
            continue
        chain_length = 1 + dfs_chain_length(acyl_atom, set([carbonyl_idx]))
        if chain_length > candidate_chain_length:
            candidate_chain_length = chain_length
    
    if candidate_chain_length == 0:
        return False, "Could not extract an acyl chain from the thioester bond"
    
    # 6. Check that the acyl chain length (including the carbonyl carbon) is in the range 6-12.
    if candidate_chain_length < 6 or candidate_chain_length > 12:
        return False, (f"Acyl chain length is {candidate_chain_length} carbons, "
                       "not in medium-chain range (6-12)")
    
    return True, (f"Matches acyl-CoA pattern with medium chain length of {candidate_chain_length} carbons "
                  "and net charge -4")

# Example usage:
if __name__ == "__main__":
    # Try one example SMILES string (cis-dec-3-enoyl-CoA(4-))
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)C/C=C\\CCCCCC)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O"
    result, reason = is_medium_chain_fatty_acyl_CoA_4__(test_smiles)
    print(result, reason)