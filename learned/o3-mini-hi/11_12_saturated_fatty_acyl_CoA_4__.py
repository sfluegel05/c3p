"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
CHEBI:84947
This classifier verifies that the molecule is a fatty acyl-CoA(4-) 
and that its fatty acyl portion – numbered linearly starting with the thioester
carbonyl as C1 – has a single (saturated) bond between C11 and C12.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-).

    The molecule must contain a fatty acyl-CoA(4-) component. We search for a short
    SMARTS fragment of CoA as a proxy and a thioester linkage ([CX3](=O)[SX2]).
    Then we locate the fatty acyl chain by starting at the carbonyl of the thioester:
    ignoring the S neighbor, we follow a single (linear) chain of connected carbon atoms.
    (If branching is encountered, the chain is not considered linear and we abort.)
    Numbering the chain such that the thioester carbonyl is C1, we require that the
    bond between C11 and C12 (i.e. between chain indices 10 and 11 in a 0-based list)
    is a single (saturated) bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): Tuple where the first element is True if the molecule meets the criteria,
                     and the second element is a reason message.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Search for a fragment indicative of a CoA moiety.
    # (This proxy is not perfect but helps ensure we are working with a fatty acyl-CoA(4-).)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found – not a fatty acyl-CoA(4-) molecule"
    
    # Look for a thioester linkage: a carbonyl directly attached to a sulfur.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester group not found – molecule may not be fatty acyl-CoA(4-)"
    
    # We assume the first thioester match is from the fatty acyl chain.
    carbonyl_idx, S_idx = thioester_matches[0][0], thioester_matches[0][1]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Find the fatty acyl chain starting from the carbonyl.
    # The carbonyl has two neighbors: one is S (we ignore) and the other (a carbon) should
    # begin the acyl chain. (The carbonyl itself is designated as C1.)
    acyl_start = None
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetIdx() != S_idx and nbr.GetAtomicNum() == 6:
            acyl_start = nbr
            break
    if acyl_start is None:
        return False, "Fatty acyl chain not found at the thioester group"
    
    # Build a linear chain (unbranched) list of carbon atoms.
    # We start with the carbonyl atom as C1 and the acyl_start as C2.
    chain = [carbonyl_atom, acyl_start]
    
    # Function to “walk” along a linear chain.
    # At each step, from the current atom, look at carbon neighbors (excluding the previous atom)
    # and if exactly one is found, that is our next chain atom.
    while True:
        current = chain[-1]
        prev = chain[-2]
        # Get all neighboring carbons excluding the one we came from.
        next_candidates = [nbr for nbr in current.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev.GetIdx()]
        if len(next_candidates) == 1:
            chain.append(next_candidates[0])
        else:
            # Either there is a branch (or end of chain) so exit the linear walk.
            break

    # Now, by convention, the chain is numbered such that:
    # chain[0] is C1 (the carbonyl), chain[1] is the first chain carbon, etc.
    # We need at least 12 carbons (so chain length >=12) to even have an 11-12 bond.
    if len(chain) < 12:
        return False, f"Fatty acyl chain is too short ({len(chain)} carbons) to have an 11-12 bond"
    
    # For numbering clarity, let C11 = chain[10] and C12 = chain[11] (0-indexed).
    atom_11 = chain[10]
    atom_12 = chain[11]
    bond_11_12 = mol.GetBondBetweenAtoms(atom_11.GetIdx(), atom_12.GetIdx())
    if bond_11_12 is None:
        return False, "11-12 bond not found in the fatty acyl chain"
    
    # Check that the 11-12 bond is a single (saturated) bond.
    if bond_11_12.GetBondType() == Chem.BondType.SINGLE:
        return True, "The fatty acyl chain has a saturated (single) bond between carbon 11 and 12"
    else:
        return False, "The 11-12 bond is not a single (saturated) bond"

# Example usage: testing with one of the provided true positive SMILES.
if __name__ == "__main__":
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_11_12_saturated_fatty_acyl_CoA_4__(test_smiles)
    print(result, ":", reason)