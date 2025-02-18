"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine
Definition:
  A 3-sn-glycerophosphoserine compound having acyl substituents at the 1- and 2-hydroxy positions.
Heuristic:
  1. Parse the SMILES and ensure at least one phosphorus atom is present.
  2. Confirm the presence of a serine fragment (using a simplified SMARTS: C(N)C(=O)O).
  3. For every ester bond (a single bond between a carbonyl carbon and an oxygen) in which:
       a. The carbonyl carbon has a double-bonded oxygen (confirming an ester carbonyl),
       b. The linking oxygen (which bridges the fatty acyl chain and the glycerol backbone)
          is connected via a short path (up to 5 bonds) to a phosphorus atom,
     follow from the carbonyl carbon the acyl chain branch (the neighbor opposite the linking oxygen)
     and count only contiguous carbon atoms.
  4. Accept the ester as valid if the acyl chain (including the carbonyl carbon) has a chain length of at least 6.
  5. Only if exactly 2 valid acyl chains are found is the molecule classified as a 3-sn-phosphatidyl-L-serine.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    Criteria:
     - The molecule must contain at least one phosphorus atom.
     - A serine fragment matching 'C(N)C(=O)O' must be present.
     - Exactly 2 ester bonds must be found in which:
         * The ester consists of a carbonyl carbon (which has a double-bonded oxygen)
           connected to a linking oxygen via a single bond.
         * The linking oxygen is connected (via a short path of at most 5 bonds) to a phosphorus atom.
         * Following the carbonyl carbon in the acyl branch (the bond not going to the linking oxygen)
           yields a contiguous chain of carbons (chain length including the carbonyl ≥ 6).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3-sn-phosphatidyl-L-serine, else False.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure at least one phosphorus (P atomic number 15) exists.
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Missing phosphorus (P) required for the phosphoglycerol head-group"
    
    # Look for a serine fragment using a simple SMARTS pattern.
    serine_smarts = "C(N)C(=O)O"
    serine_pat = Chem.MolFromSmarts(serine_smarts)
    if serine_pat is None:
        return False, "Internal error processing serine pattern"
    if not mol.HasSubstructMatch(serine_pat):
        return False, "Serine fragment (C(N)C(=O)O) not found in molecule"
    
    # Helper: Perform depth-first search to count contiguous carbons in the acyl chain.
    def dfs_chain_length(start_idx, banned):
        max_length = 0
        stack = [(start_idx, 0, set())]  # (current atom index, current chain length, visited set)
        while stack:
            cur_idx, cur_length, visited = stack.pop()
            if cur_idx in visited:
                continue
            visited = visited | {cur_idx}
            cur_length += 1
            if cur_length > max_length:
                max_length = cur_length
            atom = mol.GetAtomWithIdx(cur_idx)
            # Traverse only to carbon neighbors that are not banned.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in banned:
                    stack.append((nbr.GetIdx(), cur_length, visited))
        return max_length

    # Helper: Check if a given oxygen is connected to a phosphorus atom via a path with at most max_bonds.
    def oxygen_connected_to_phosphorus(o_atom, max_bonds=5):
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 15:  # phosphorus
                try:
                    # GetShortestPath returns a tuple of atom indices
                    path = Chem.GetShortestPath(mol, o_atom.GetIdx(), atom.GetIdx())
                    if path and (len(path) - 1) <= max_bonds:
                        return True
                except Exception:
                    continue
        return False

    valid_acyl_count = 0
    counted_esters = set()  # To avoid double-counting same ester bond.
    
    # Iterate over all bonds to detect candidate ester bonds.
    for bond in mol.GetBonds():
        # Consider only single bonds (the linking bond in an ester is single).
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        
        # Identify a candidate: one atom must be a carbon (potential carbonyl) and the other oxygen.
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8:
            carbonyl = a1
            linking_oxygen = a2
        elif a2.GetAtomicNum() == 6 and a1.GetAtomicNum() == 8:
            carbonyl = a2
            linking_oxygen = a1
        else:
            continue
        
        # Confirm that the carbonyl carbon has a double-bonded oxygen.
        found_double = False
        for nbr in carbonyl.GetNeighbors():
            b = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and b is not None and b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                found_double = True
                break
        if not found_double:
            continue
        
        # Verify that the linking oxygen is connected to phosphorus with a path of ≤5 bonds.
        if not oxygen_connected_to_phosphorus(linking_oxygen, max_bonds=5):
            continue
        
        # From the carbonyl carbon, select the neighbor that is not the linking oxygen and not the double-bonded oxygen.
        acyl_start = None
        for nbr in carbonyl.GetNeighbors():
            if nbr.GetIdx() == linking_oxygen.GetIdx():
                continue
            b = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
            # Skip the double-bonded oxygen.
            if nbr.GetAtomicNum() == 8 and b is not None and b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_start = nbr
                break
        if acyl_start is None:
            continue
        
        # Use a unique identifier for the ester bond to avoid double-counting.
        bond_id = tuple(sorted((carbonyl.GetIdx(), linking_oxygen.GetIdx())))
        if bond_id in counted_esters:
            continue
        counted_esters.add(bond_id)
        
        # Count the acyl chain length: include the carbonyl carbon plus contiguous carbons from the acyl branch.
        chain_length = 1 + dfs_chain_length(acyl_start.GetIdx(), banned={carbonyl.GetIdx()})
        if chain_length >= 6:
            valid_acyl_count += 1

    if valid_acyl_count != 2:
        return False, f"Expected 2 acyl groups with chain length ≥6 in the ester bonds anchored to the phosphate backbone; found {valid_acyl_count}"
    
    return True, "Molecule contains a phosphoserine head-group with exactly 2 acyl groups anchored to the phosphoglycerol backbone"

# Example usage (for testing – these may be commented out in production):
if __name__ == "__main__":
    examples = [
        ("1,2-distearoyl-sn-glycero-3-phosphoserine", "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCCCCCCCCCCCC"),
        ("1-oleoyl-2-arachidonoyl-sn-glycero-3-phospho-L-serine", "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC")
    ]
    for name, smi in examples:
        result, reason = is_3_sn_phosphatidyl_L_serine(smi)
        print(name, "->", result, reason)