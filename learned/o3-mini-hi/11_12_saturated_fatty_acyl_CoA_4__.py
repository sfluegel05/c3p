"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
CHEBI:84947

This classifier verifies that the molecule is a fatty acyl-CoA(4-)
and that its fatty acyl portion – numbered linearly starting with the thioester
carbonyl as C1 – has a saturated (single) bond between C11 and C12.
It uses a DFS strategy to determine the longest carbon chain emanating from
the thioester carbonyl (ignoring the sulfur) so that branching (e.g. methyl branches)
does not lead prematurely to aborting the acyl chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-).

    The molecule must contain a fatty acyl-CoA(4-) component. First we search
    for a proxy SMARTS for CoA and a thioester linkage. Then, starting from the
    thioester carbonyl we find all possible linear acyl chains (by a DFS that only
    follows carbon neighbors) and choose the longest chain. Numbering the chain
    such that the thioester carbonyl is C1, we then require that the bond between
    C11 and C12 (i.e. chain indices 10 and 11 in zero-based indexing) is a single bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): Tuple where the first element is True if the molecule meets the criteria,
                     and the second element a reason message.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a proxy pattern for a CoA fragment.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found – not a fatty acyl-CoA(4-) molecule"
    
    # Look for a thioester linkage: a carbonyl ([CX3](=O)) directly attached to a sulfur ([SX2])
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester group not found – molecule may not be fatty acyl-CoA(4-)"
    
    # Assume the first thioester match corresponds to the fatty acyl chain.
    # thioester_matches returns tuples of atom indices; in our SMARTS, the first is carbonyl, second is S.
    carbonyl_idx, sulfur_idx = thioester_matches[0][0], thioester_matches[0][1]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # The carbonyl should have two neighbors. One is S (thioester sulfur) which we ignore.
    # The other, if it is carbon, starts the fatty acyl chain.
    acyl_start = None
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetIdx() != sulfur_idx and nbr.GetAtomicNum() == 6:
            acyl_start = nbr
            break
    if acyl_start is None:
        return False, "Fatty acyl chain not found at the thioester group"
    
    # To follow the chain we perform a DFS from the acyl_start.
    # We only traverse atoms that are carbons and do not revisit atoms.
    def dfs(current_atom, visited):
        # Each path is a list of atom objects - starting with current_atom.
        paths = [[current_atom]]
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:  # only follow carbon atoms
                continue
            if nbr.GetIdx() in visited:
                continue
            # Continue the DFS from this neighbor.
            new_visited = visited | {nbr.GetIdx()}
            sub_paths = dfs(nbr, new_visited)
            for path in sub_paths:
                paths.append([current_atom] + path)
        return paths

    # We begin the DFS from acyl_start with a visited set containing its index.
    all_paths = dfs(acyl_start, {acyl_start.GetIdx()})
    # From these possible paths, choose the one with maximal length.
    # (If there are several of equal length, just take the first.)
    if not all_paths:
        return False, "No acyl chain path found"
    longest_path = max(all_paths, key=len)

    # Our fatty acyl chain is defined as the thioester carbonyl followed by the longest path.
    # Numbering: chain[0] is thioester carbonyl (C1), chain[1] is first carbon (C2), etc.
    acyl_chain = [carbonyl_atom] + longest_path
    
    # A valid fatty acyl chain must have at least 12 carbons to contain an 11-12 bond.
    if len(acyl_chain) < 12:
        return False, f"Fatty acyl chain is too short ({len(acyl_chain)} atoms including the carbonyl) to have an 11-12 bond"
    
    # By numbering, C11 is acyl_chain[10] and C12 is acyl_chain[11].
    atom_11 = acyl_chain[10]
    atom_12 = acyl_chain[11]
    bond_11_12 = mol.GetBondBetweenAtoms(atom_11.GetIdx(), atom_12.GetIdx())
    if bond_11_12 is None:
        return False, "11-12 bond not found in the fatty acyl chain"
    
    # Check that the 11-12 bond is a single (saturated) bond. (Double or triple bonds indicate unsaturation.)
    if bond_11_12.GetBondType() == Chem.BondType.SINGLE:
        return True, "The fatty acyl chain has a saturated (single) bond between carbon 11 and 12"
    else:
        return False, "The 11-12 bond is not a single (saturated) bond"

# Example usage: testing with one of the provided true positive SMILES.
if __name__ == "__main__":
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_11_12_saturated_fatty_acyl_CoA_4__(test_smiles)
    print(result, ":", reason)