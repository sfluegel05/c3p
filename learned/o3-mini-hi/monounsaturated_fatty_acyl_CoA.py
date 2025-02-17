"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: Monounsaturated fatty acyl-CoA
Definition: Any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.
Improved strategy:
  - Identify the thioester linkage [CX3](=O)[SX2].
  - From the carbonyl carbon, choose the substituent that is a carbon (this becomes the start of the fatty acyl chain).
  - Perform a DFS to collect only aliphatic (non-aromatic, non-ring) carbon atoms – a heuristic to isolate the fatty acyl chain.
  - Check that the chain has at least 6 carbons.
  - Count only bonds between chain atoms that are carbon–carbon double bonds.
  - Ensure exactly one such C=C bond is found.
  - Also verify that the molecule contains a CoA moiety by matching the adenine fragment.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA.
    That is, the molecule must contain:
      (a) a thioester bond in the form of [CX3](=O)[SX2],
      (b) a CoA moiety (detected via the adenine substructure),
      (c) a fatty acyl chain (the chain attached to the carbonyl carbon) that is aliphatic,
          at least 6 carbons long, and contains exactly one carbon-carbon double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a monounsaturated fatty acyl-CoA, False otherwise.
        str: A reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the thioester group: [C](=O)[S]
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found (lack of fatty acyl linkage)"
    
    # Assume the first thioester match corresponds to the linkage.
    # In the match, atom index 0 is the carbonyl carbon.
    thioester_match = thioester_matches[0]
    carbonyl_idx = thioester_match[0]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Identify the fatty acyl chain substituent:
    # The carbonyl carbon is connected to:
    #  - an oxygen (double bond) and
    #  - a sulfur (leading to CoA)
    # So the remaining neighbor (if any) should be the first carbon in the acyl chain.
    alpha_atom = None
    for neighbor in carbonyl_atom.GetNeighbors():
        # Exclude oxygen and sulfur.
        if neighbor.GetAtomicNum() == 6:
            alpha_atom = neighbor
            break
    if alpha_atom is None:
        return False, "No fatty acyl chain detected from thioester (missing R-group)"
    
    # Perform a DFS to extract the fatty acyl chain.
    # We restrict to non-aromatic, acyclic carbons (heuristic for an aliphatic chain)
    chain_atom_indices = set()
    stack = [carbonyl_atom]  # include the carbonyl carbon as part of the chain
    while stack:
        atom = stack.pop()
        idx = atom.GetIdx()
        # Only consider carbon atoms
        if atom.GetAtomicNum() != 6:
            continue
        # Exclude aromatic atoms or atoms in rings (to avoid picking up conjugated/aromatic groups)
        if atom.GetIsAromatic() or atom.IsInRing():
            continue
        # Mark as part of fatty acyl chain
        if idx in chain_atom_indices:
            continue
        chain_atom_indices.add(idx)
        # Traverse neighbors but do not cross into heteroatoms (like O, S, N) or CoA regions.
        for nbr in atom.GetNeighbors():
            # Only follow if neighbor is carbon.
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
                if nbr.GetIdx() not in chain_atom_indices:
                    stack.append(nbr)
    
    # Minimum chain length check (typically a fatty acyl chain is long)
    if len(chain_atom_indices) < 6:
        return False, f"Extracted fatty acyl chain is too short ({len(chain_atom_indices)} carbons)"
    
    # Count carbon-carbon double bonds within the extracted fatty chain.
    double_bond_count = 0
    for bond in mol.GetBonds():
        # Check if both atoms belong to the chain.
        idx1 = bond.GetBeginAtom().GetIdx()
        idx2 = bond.GetEndAtom().GetIdx()
        if idx1 in chain_atom_indices and idx2 in chain_atom_indices:
            # Only count bonds where both atoms are carbons.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1

    if double_bond_count != 1:
        return False, f"Fatty acyl chain has {double_bond_count} carbon-carbon double bond(s); exactly one is required"
    
    # Verify that the molecule contains a CoA moiety.
    # We look for the adenine substructure which appears as: N1C=NC2=C1N=CN=C2
    coa_pattern = Chem.MolFromSmarts("N1C=NC2=C1N=CN=C2")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety (adenine fragment) not found"
    
    return True, "Contains an aliphatic fatty acyl chain with exactly one C=C double bond and a CoA moiety"

# Example usage: testing with one of the provided SMILES (oleoyl-CoA)
if __name__ == "__main__":
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N"
    result, reason = is_monounsaturated_fatty_acyl_CoA(test_smiles)
    print(result, reason)