"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
CHEBI:84947 is mentioned in the definition.
This classifier verifies that the molecule is a fatty acyl-CoA(4-) 
and that the fatty acyl portion (the longest carbon chain starting from the thioester carbonyl) 
has a saturated (i.e. single) bond between C11 and C12 (numbering from the carbonyl carbon as C1).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-).

    A fatty acyl-CoA(4-) contains a fatty acyl chain attached via a thioester 
    (i.e. a [CX3](=O)[SX2] functional group) to a CoA moiety. 
    The 11-12-saturated fatty acyl-CoA(4-) class is defined as any such molecule in which
    the bond connecting the 11th and 12th carbon atoms (counting from the carbonyl carbon as C1)
    in the fatty acyl chain is saturated (a single bond).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): (True, reason) if molecule matches the class,
                     (False, reason) if not.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Verify presence of CoA moiety.
    # We search for a fragment of the CoA tail; note that this is an imperfect proxy.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found - not a fatty acyl-CoA(4-) molecule"
    
    # Look for thioester linkage: a carbonyl directly attached to a sulfur.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester group not found - molecule may not be fatty acyl-CoA(4-)"
    
    # We assume the first thioester match is from the fatty acyl chain.
    carbonyl_idx, S_idx = thioester_matches[0][0], thioester_matches[0][1]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # From the carbonyl, one neighbor is S (which we ignore). The other carbon neighbor
    # should begin the fatty acyl chain.
    acyl_start = None
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetIdx() != S_idx and nbr.GetAtomicNum() == 6:
            acyl_start = nbr
            break
    if acyl_start is None:
        return False, "Fatty acyl chain not found at the thioester group"
    
    # Define a helper function that, given a starting atom, recursively obtains
    # the longest simple (non-cyclic) path made only of carbon atoms.
    def get_longest_chain_path(atom, visited):
        # Mark current atom as visited
        visited = visited | {atom.GetIdx()}
        path = [atom]
        longest_suffix = []
        for nbr in atom.GetNeighbors():
            # Only follow carbon atoms not yet visited.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                candidate = get_longest_chain_path(nbr, visited)
                if len(candidate) > len(longest_suffix):
                    longest_suffix = candidate
        return path + longest_suffix

    # We now find the longest acyl chain starting from the carbonyl.
    # According to convention, carbonyl carbon (from the thioester) is C1.
    # From carbonyl take the branch from the neighbor that is not the sulfur.
    # But to cover branching in the chain, we consider paths starting from carbonyl.
    longest_chain = get_longest_chain_path(carbonyl_atom, set())
    
    # Ensure the chain represents a fatty acyl group by requiring at least 12 carbons;
    # note that the chain includes the carbonyl carbon as C1.
    if len(longest_chain) < 12:
        return False, f"Fatty acyl chain is too short (only {len(longest_chain)} carbons) to have an 11-12 bond"
    
    # Identify the atoms corresponding to positions C11 and C12.
    # (0-indexed: chain index 10 is C11, and index 11 is C12)
    atom_11 = longest_chain[10]
    atom_12 = longest_chain[11]
    bond_11_12 = mol.GetBondBetweenAtoms(atom_11.GetIdx(), atom_12.GetIdx())
    
    if bond_11_12 is None:
        return False, "11-12 bond not found in the fatty acyl chain"
    
    # Check that the bond is a single (saturated) bond.
    if bond_11_12.GetBondType() == Chem.BondType.SINGLE:
        return True, "The fatty acyl chain has a saturated (single) bond between carbon 11 and 12"
    else:
        return False, "The 11-12 bond is unsaturated (not a single bond), which does not meet the class criteria"

# Example usage: testing with a provided true positive
if __name__ == "__main__":
    # Example test: (13Z)-3-oxodocosenoyl-CoA(4-) from the provided list
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_11_12_saturated_fatty_acyl_CoA_4__(test_smiles)
    print(result, ":", reason)