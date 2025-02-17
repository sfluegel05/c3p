"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: Monounsaturated fatty acyl-CoA

Definition: Any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon–carbon double bond.

Strategy:
  - Look for the thioester linkage [CX3](=O)[SX2] linking the fatty acyl chain to the CoA.
  - From the carbonyl carbon, identify the carbon neighbor (the α‐carbon) from which the fatty acyl chain extends.
  - Traverse the fatty acyl chain (following non-aromatic, non-ring carbon atoms) and count its carbons.
  - Count the number of carbon–carbon double bonds in that chain (exactly one is required).
  - Verify the presence of the CoA moiety by searching for an adenine-like substructure. 
    To improve the matching, stereochemistry will be removed from the molecule and multiple adenine SMARTS patterns will be used.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA.
    It requires the molecule to have:
      (a) a thioester linkage [CX3](=O)[SX2],
      (b) an aliphatic fatty acyl chain (starting at the alpha carbon off the carbonyl) that is at
          least 6 carbons long and contains exactly one carbon–carbon double bond,
      (c) a CoA moiety which is identified by the presence of an adenine-like substructure.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if all criteria are met, otherwise False.
        str: A reason string explaining the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Remove stereochemistry to ease matching of adenine substructure.
    Chem.RemoveStereochemistry(mol)
    
    # --- 1. Locate the thioester linkage: [CX3](=O)[SX2]
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found (fatty acyl linkage missing)"
    
    # Choose the first match (assuming it represents the fatty acyl thioester).
    thioester_match = thioester_matches[0]
    # The first atom in the match is the carbonyl carbon.
    carbonyl_idx = thioester_match[0]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # --- 2. Identify the α-carbon in the acyl chain. 
    # This is the neighbor to the carbonyl carbon that is a carbon (not sulfur).
    alpha_atom = None
    for neighbor in carbonyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            alpha_atom = neighbor
            break
    if alpha_atom is None:
        return False, "No fatty acyl chain detected from thioester linkage (missing R-group)"
    
    # --- 3. Traverse the acyl chain via depth-first search (DFS) starting from the α-carbon.
    # Only follow carbon atoms that are not aromatic and not part of any ring.
    chain_atom_indices = set()
    stack = [alpha_atom]
    while stack:
        atom = stack.pop()
        idx = atom.GetIdx()
        if atom.GetAtomicNum() != 6:
            continue
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
        # Count the bond if both atoms are in the chain and the bond is a double bond.
        if a1.GetIdx() in chain_atom_indices and a2.GetIdx() in chain_atom_indices:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1
    if double_bond_count != 1:
        return False, f"Fatty acyl chain has {double_bond_count} C=C double bond(s); exactly one is required"
    
    # --- 5. Verify the presence of the CoA moiety by detecting an adenine-like fragment.
    # Use multiple SMARTS patterns to capture different representations of adenine.
    adenine_smarts_list = [
        "N1C=NC2=C1N=CN=C2",   # common representation
        "n1cnc2ncnc(n2)n1"      # alternative pattern ignoring stereochem details
    ]
    adenine_found = False
    for smarts in adenine_smarts_list:
        adenine_pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(adenine_pattern):
            adenine_found = True
            break
    if not adenine_found:
        return False, "CoA moiety (adenine fragment) not found"
    
    return True, "Contains fatty acyl chain with exactly one C=C double bond and a CoA moiety"

# Example usage: testing the function with one of the provided SMILES strings.
if __name__ == "__main__":
    # Test with (E)-hexadec-2-enoyl-CoA.
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_monounsaturated_fatty_acyl_CoA(test_smiles)
    print(result, reason)