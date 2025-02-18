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
       a. The carbonyl carbon has a double-bonded oxygen (confirming ester carbonyl),
       b. The oxygen (the ester oxygen linking the acyl chain) is “anchored” via a short path 
          (shortest path length ≤ 4 bonds) to at least one phosphorus atom,
     follow from the carbonyl carbon the acyl chain branch (the neighbor opposite the ester oxygen)
     and count only contiguous carbon atoms.
  4. Accept the ester as valid if the acyl chain has a chain length (including the carbonyl carbon)
     of at least 6.
  5. Only if exactly 2 valid acyl chains are found is the molecule classified as a 3-sn-phosphatidyl-L-serine.
"""
from rdkit import Chem
from collections import deque

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    The criteria are:
     - The molecule must contain at least one phosphorus atom.
     - A serine fragment matching 'C(N)C(=O)O' (ignoring stereochemistry) must be present.
     - There must be exactly 2 ester bonds where:
         * The ester motif involves a carbonyl carbon (with a double-bonded oxygen) 
           and a linking oxygen;
         * The linking oxygen is connected (via a path of at most 4 bonds) to a phosphorus atom;
         * Following the carbonyl carbon toward the acyl chain (the bond not going to the linking oxygen)
           yields a chain made only of carbons having a length of at least 6 atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3-sn-phosphatidyl-L-serine, else False.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that at least one phosphorus exists (P atomic number 15)
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Missing phosphorus (P) required for the phosphoglycerol head-group"
    
    # Check for serine fragment using a simple SMARTS pattern
    serine_smarts = "C(N)C(=O)O"
    serine_pat = Chem.MolFromSmarts(serine_smarts)
    if serine_pat is None:
        return False, "Internal error processing serine pattern"
    if not mol.HasSubstructMatch(serine_pat):
        return False, "Serine fragment (C(N)C(=O)O) not found in molecule"
    
    # Helper: Given a starting carbon (by its index), perform a depth-first search counting contiguous carbons.
    def dfs_chain_length(start_idx, banned):
        max_length = 0
        stack = [(start_idx, 0, set())]  # (current_atom_idx, current_length, visited)
        while stack:
            cur_idx, cur_length, visited = stack.pop()
            if cur_idx in visited:
                continue
            visited = visited | {cur_idx}
            cur_length += 1
            if cur_length > max_length:
                max_length = cur_length
            atom = mol.GetAtomWithIdx(cur_idx)
            # traverse only to carbon neighbors not banned (avoid going back into backbone)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in banned:
                    stack.append((nbr.GetIdx(), cur_length, visited))
        return max_length
    
    # Helper: Check if a given atom (typically the ester oxygen) is connected to a phosphorus atom
    # via a shortest path of length (number of bonds) less than or equal to max_bonds.
    def oxygen_connected_to_phosphorus(o_atom, max_bonds=4):
        # iterate over all phosphorus atoms and compute shortest path lengths
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 15:
                path = Chem.GetShortestPath(mol, o_atom.GetIdx(), atom.GetIdx())
                if path and (len(path)-1) <= max_bonds:
                    return True
        return False
    
    valid_acyl_count = 0
    counted_esters = set()  # avoid double counting the same ester bond
    
    # Iterate over all bonds to detect candidate ester bonds.
    for bond in mol.GetBonds():
        # Only consider single bonds (the ester link between carbonyl carbon and O)
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        
        # Identify candidate where one atom is carbon and the other is oxygen.
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8:
            carbonyl = a1
            estero = a2  # candidate linking oxygen
        elif a2.GetAtomicNum() == 6 and a1.GetAtomicNum() == 8:
            carbonyl = a2
            estero = a1
        else:
            continue
        
        # Check that the carbonyl carbon has a double-bonded oxygen (carbonyl O)
        found_dbl = False
        for nbr in carbonyl.GetNeighbors():
            b = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and b is not None and b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                found_dbl = True
                break
        if not found_dbl:
            continue
        
        # Verify that the linking oxygen (estero) is anchored to the phosphate
        if not oxygen_connected_to_phosphorus(estero, max_bonds=4):
            continue
        
        # Identify acyl branch from carbonyl: choose the neighbor that is NOT the linking oxygen
        # and not the double-bonded O.
        acyl_start = None
        for nbr in carbonyl.GetNeighbors():
            # skip the ester oxygen
            if nbr.GetIdx() == estero.GetIdx():
                continue
            # skip the double-bonded oxygen
            b = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and b is not None and b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_start = nbr
                break
        if acyl_start is None:
            continue
        
        # Use unique identifier for the ester bond to avoid double-counting.
        bond_id = tuple(sorted((carbonyl.GetIdx(), estero.GetIdx())))
        if bond_id in counted_esters:
            continue
        counted_esters.add(bond_id)
        
        # Count the carbon chain length. We count the carbonyl carbon (1)
        # plus the contiguous chain in the acyl branch, avoiding going back into the carbonyl.
        chain_length = 1 + dfs_chain_length(acyl_start.GetIdx(), banned={carbonyl.GetIdx()})
        if chain_length >= 6:
            valid_acyl_count += 1

    if valid_acyl_count != 2:
        return False, f"Expected 2 acyl groups with chain length ≥6 in the ester bonds anchored to the phosphate backbone; found {valid_acyl_count}"
    
    return True, "Molecule contains a phosphoserine head-group and exactly 2 acyl groups (chain length ≥6) anchored to a phosphoglycerol backbone"

# Example usage (for testing; these can be commented out):
if __name__ == "__main__":
    examples = [
        ("1,2-distearoyl-sn-glycero-3-phosphoserine", "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCCCCCCCCCCCC"),
        ("1-oleoyl-2-arachidonoyl-sn-glycero-3-phospho-L-serine", "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC")
    ]
    for name, smi in examples:
        result, reason = is_3_sn_phosphatidyl_L_serine(smi)
        print(name, "->", result, reason)