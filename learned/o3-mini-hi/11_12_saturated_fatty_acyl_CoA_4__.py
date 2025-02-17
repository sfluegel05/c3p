"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
CHEBI:84947

This classifier checks that the molecule contains a fatty acyl-CoA(4-) unit.
It first locates a CoA fragment (using a proxy SMARTS) and a thioester group.
Then it “isolates” the fatty acyl chain by excluding atoms that belong to the CoA moiety
(as determined by the SMARTS match) and following only carbon atoms from the carbonyl.
It then computes the longest simple (no revisiting) carbon-only chain starting from the thioester
carbonyl’s carbon neighbor (ignoring the sulfur connection).
Finally, it numbers the chain with the carbonyl as C1 and checks that the bond
between C11 and C12 (atoms at positions 10 and 11) is a single bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-).
    
    The classifier requires:
      • A CoA fragment (using a proxy SMARTS: "SCCNC(=O)CCNC(=O)")
      • A thioester linkage ([CX3](=O)[SX2])
      • The fatty acyl chain (taken as the longest carbon‐only path emanating from the
        thioester carbonyl, while excluding atoms belonging to the CoA portion)
      • That the bond between positions C11 and C12 (numbered with the carbonyl as C1)
        is a single (saturated) bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple in which the first element is True if the molecule conforms
                     to the class; otherwise False. The second element is a message giving
                     the reason for the classification.
    """
    # Get the rdkit molecule:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a CoA fragment by a proxy SMARTS.
    # Note: the pattern here is a proxy for part of CoA.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "CoA moiety not found – not a fatty acyl-CoA(4-) molecule"
    # We’ll exclude atoms in any (first) CoA match from our acyl chain search.
    coa_indices = set(coa_matches[0])
    
    # Look for a thioester linkage: a carbonyl ([CX3](=O)) directly attached to a sulfur ([SX2]).
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester group not found – molecule may not be fatty acyl-CoA(4-)"
    # We assume the first match is the fatty acyl thioester.
    carbonyl_idx, sulfur_idx = thioester_matches[0][0], thioester_matches[0][1]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Identify the neighbor of the carbonyl atom that will be the start of the acyl chain.
    # We skip the sulfur and any atom belonging to the CoA pattern.
    acyl_start = None
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetIdx() == sulfur_idx:
            continue
        if nbr.GetIdx() in coa_indices:
            continue
        if nbr.GetAtomicNum() == 6:
            acyl_start = nbr
            break
    if acyl_start is None:
        return False, "Fatty acyl chain not found at the thioester group"
    
    # To “isolate” the fatty acyl chain we will perform a DFS starting from the acyl_start atom.
    # We follow only carbon atoms (atomic number 6) and we exclude atoms that belong to the CoA moiety.
    # Our goal is to find the longest simple (no repeated atoms) path.
    def longest_chain(atom, visited):
        best_path = [atom]
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in visited:
                continue
            if nbr.GetIdx() in coa_indices:
                continue
            new_visited = visited | {nbr.GetIdx()}
            candidate = longest_chain(nbr, new_visited)
            # candidate is path starting at nbr; add current atom at front.
            candidate_path = [atom] + candidate
            if len(candidate_path) > len(best_path):
                best_path = candidate_path
        return best_path

    # Our complete fatty acyl chain is defined as:
    # [carbonyl_atom] followed by the longest chain starting from acyl_start.
    chain = [carbonyl_atom] + longest_chain(acyl_start, {acyl_start.GetIdx()})
    
    # The fatty acyl chain numbering: C1 is the thioester carbonyl; C2, C3,... follow from the chain.
    # To have a bond between C11 and C12, we need at least 12 atoms:
    if len(chain) < 12:
        return False, f"Fatty acyl chain is too short ({len(chain)} atoms including the carbonyl) to have an 11-12 bond"
    
    # Check the bond type between C11 (chain index 10) and C12 (chain index 11).
    atom_11 = chain[10]
    atom_12 = chain[11]
    bond = mol.GetBondBetweenAtoms(atom_11.GetIdx(), atom_12.GetIdx())
    if bond is None:
        return False, "No bond found between potential C11 and C12 atoms"
    if bond.GetBondType() == Chem.BondType.SINGLE:
        return True, "The fatty acyl chain has a saturated (single) bond between carbon 11 and 12"
    else:
        return False, f"The bond between carbon 11 and 12 is {bond.GetBondType()}, not a single bond"

# Example usage:
if __name__ == "__main__":
    # Test with a true positive:
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_11_12_saturated_fatty_acyl_CoA_4__(test_smiles)
    print(result, ":", reason)