"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: Monounsaturated fatty acyl-CoA
Definition: Any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon–carbon double bond.
Improved strategy:
  - Identify the thioester group [CX3](=O)[SX2].
  - From the carbonyl carbon, choose the neighbor that is a carbon (this becomes the starting α-carbon for the fatty acyl chain).
  - Perform a DFS starting from the α-carbon to collect only non‐aromatic, acyclic carbons (heuristic to isolate an aliphatic fatty acyl chain).
  - Check that the chain has at least 6 carbon atoms.
  - Count the carbon–carbon double bonds (only between atoms in the fatty chain) and require exactly one.
  - Verify that the molecule contains a CoA moiety by matching the adenine fragment using the SMARTS pattern "n1cnc2ncnc12".
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA.
    Specifically, the molecule must contain:
      (a) a thioester bond in the form of [CX3](=O)[SX2],
      (b) a fatty acyl chain (connected to the carbonyl carbon via an α-carbon) that is aliphatic,
          at least 6 carbons long, and contains exactly one carbon–carbon double bond,
      (c) a CoA moiety as identified by the adenine substructure.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule satisfies all the criteria, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- 1. Identify thioester group: [CX3](=O)[SX2]
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found (fatty acyl linkage missing)"
    
    # Assume first match corresponds to fatty acyl thioester linkage.
    # In the match: atom index 0 is the carbonyl carbon.
    thioester_match = thioester_matches[0]
    carbonyl_idx = thioester_match[0]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # --- 2. Identify the fatty acyl chain connection
    # The carbonyl is attached to an oxygen (double bond O) and a sulfur (attached to CoA).
    # So the remaining neighbor (a carbon) is the start (alpha carbon) of the fatty acyl chain.
    alpha_atom = None
    for neighbor in carbonyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            alpha_atom = neighbor
            break
    if alpha_atom is None:
        return False, "No fatty acyl chain detected from thioester linkage (missing R-group)"
    
    # --- 3. Extract the fatty acyl chain using DFS from the alpha carbon.
    # We only follow non-aromatic, non-ring carbon atoms.
    chain_atom_indices = set()
    stack = [alpha_atom]
    while stack:
        atom = stack.pop()
        idx = atom.GetIdx()
        # Only consider carbon atoms
        if atom.GetAtomicNum() != 6:
            continue
        # Exclude aromatic atoms or atoms in rings (to avoid conjugated or cyclic portions).
        if atom.GetIsAromatic() or atom.IsInRing():
            continue
        if idx in chain_atom_indices:
            continue
        chain_atom_indices.add(idx)
        # Traverse neighbors that are carbons (and meet our aliphatic criteria).
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
                if nbr.GetIdx() not in chain_atom_indices:
                    stack.append(nbr)
    
    if len(chain_atom_indices) < 6:
        return False, f"Extracted fatty acyl chain is too short ({len(chain_atom_indices)} carbons found)"
    
    # --- 4. Count the carbon–carbon double bonds within the extracted chain.
    double_bond_count = 0
    # Iterate over bonds in the molecule; consider bonds connecting two atoms in the chain.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetIdx() in chain_atom_indices and a2.GetIdx() in chain_atom_indices:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1
    if double_bond_count != 1:
        return False, f"Fatty acyl chain has {double_bond_count} carbon–carbon double bond(s); exactly one is required"
    
    # --- 5. Verify presence of the CoA moiety by looking for the adenine fragment.
    # Use a more general SMARTS pattern for adenine using aromatic atoms.
    adenine_pattern = Chem.MolFromSmarts("n1cnc2ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "CoA moiety (adenine fragment) not found"
    
    return True, "Contains an aliphatic fatty acyl chain with exactly one C=C double bond and a CoA moiety"

# Example usage: you can test with one of the provided SMILES strings.
if __name__ == "__main__":
    # Example: (E)-hexadec-2-enoyl-CoA
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_monounsaturated_fatty_acyl_CoA(test_smiles)
    print(result, reason)