"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: Monounsaturated fatty acyl-CoA

Definition: Any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon–carbon double bond.
Strategy:
  - Identify the thioester group [CX3](=O)[SX2] which links the fatty acyl chain to CoA.
  - From the carbonyl carbon, pick the attached carbon (the α-carbon of the acyl chain).
  - Extract the acyl chain via DFS, following only aliphatic (non-aromatic, non‐ring) carbon atoms.
  - Ensure the chain has at least 6 carbon atoms, and exactly one C=C double bond.
  - Verify that the structure contains a CoA moiety by searching for an adenine fragment.
    We use two SMARTS patterns ("n1cnc2ncnc12" and "N1C=NC2=C1N=CN=C2") to robustly detect the adenine portion.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA.
    It requires the molecule to have:
      (a) a thioester linkage [CX3](=O)[SX2],
      (b) an aliphatic fatty acyl chain (from the alpha carbon off the carbonyl) that is at
          least 6 carbons long and contains exactly one carbon–carbon double bond,
      (c) a CoA moiety as shown by an adenine substructure.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if all criteria are met.
        str: A reason string explaining the decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- 1. Identify thioester linkage: [CX3](=O)[SX2]
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found (fatty acyl linkage missing)"
    
    # Assume the first match corresponds to the fatty acyl thioester linkage.
    # In the match, the first atom (index 0) is the carbonyl carbon.
    thioester_match = thioester_matches[0]
    carbonyl_idx = thioester_match[0]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # --- 2. Identify the alpha (α) carbon that initiates the fatty acyl chain.
    # The carbonyl group is bound to an oxygen (as a carbonyl) and a sulfur (attached to CoA).
    # So the remaining neighbor that is a carbon is the fatty chain's starting point.
    alpha_atom = None
    for neighbor in carbonyl_atom.GetNeighbors():
        # We expect one neighbor to be sulfur (of the thioester) and one to be a carbon.
        if neighbor.GetAtomicNum() == 6:
            alpha_atom = neighbor
            break
    if alpha_atom is None:
        return False, "No fatty acyl chain detected from thioester linkage (missing R-group)"
    
    # --- 3. Extract the fatty acyl chain via DFS from the alpha carbon.
    # Only follow non-aromatic, non-ring carbon atoms.
    chain_atom_indices = set()
    stack = [alpha_atom]
    while stack:
        atom = stack.pop()
        idx = atom.GetIdx()
        if atom.GetAtomicNum() != 6:  # only carbon atoms
            continue
        # Exclude aromatic or ring atoms (to ignore conjugated or cyclic parts)
        if atom.GetIsAromatic() or atom.IsInRing():
            continue
        if idx in chain_atom_indices:
            continue
        chain_atom_indices.add(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
                if nbr.GetIdx() not in chain_atom_indices:
                    stack.append(nbr)
    
    if len(chain_atom_indices) < 6:
        return False, f"Extracted fatty acyl chain is too short ({len(chain_atom_indices)} carbons found)"
    
    # --- 4. Count the carbon–carbon double bonds within the extracted fatty acyl chain.
    double_bond_count = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetIdx() in chain_atom_indices and a2.GetIdx() in chain_atom_indices:
            # Count only actual double bonds between carbons.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1
    if double_bond_count != 1:
        return False, f"Fatty acyl chain has {double_bond_count} carbon–carbon double bond(s); exactly one is required"
    
    # --- 5. Verify presence of the CoA moiety.
    # Use two adenine-related SMARTS patterns to increase detection robustness.
    adenine_smarts_list = ["n1cnc2ncnc12", "N1C=NC2=C1N=CN=C2"]
    adenine_found = False
    for smarts in adenine_smarts_list:
        adenine_pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(adenine_pattern):
            adenine_found = True
            break
    if not adenine_found:
        return False, "CoA moiety (adenine fragment) not found"
    
    return True, "Contains aliphatic fatty acyl chain with exactly one C=C double bond and a CoA moiety"

# Example usage: test with one of the provided SMILES strings.
if __name__ == "__main__":
    # Example: (E)-hexadec-2-enoyl-CoA
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_monounsaturated_fatty_acyl_CoA(test_smiles)
    print(result, reason)