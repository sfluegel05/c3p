"""
Classifies: CHEBI:17389 2-monoglyceride
"""
#!/usr/bin/env python3
"""
Classifies: 2-monoglyceride

A 2-monoglyceride is defined as a monoglyceride in which the acyl substituent is located 
at the 2-position. That is, the glycerol backbone is modified from HOCH2–CH(OH)–CH2OH 
into HOCH2–CH(O–C(=O)R)–CH2OH.

This program uses a SMARTS substructure to identify the key fragment:
    HOCH2–CH(O–C(=O)R)–CH2OH

If a valid fragment is found and the acyl chain (R) has at least 5 contiguous aliphatic carbons,
the molecule is classified as a 2-monoglyceride.
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    
    A 2-monoglyceride contains a glycerol backbone where the acyl (ester) group is attached 
    at the central (2-) position. That is, one of the –OH groups of glycerol is replaced by an acyl ester.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to get reliable connectivity
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS query for the 2-monoglyceride motif.
    # This pattern aims to match a glycerol fragment with two terminal CH2OH groups and 
    # an esterified group at the central carbon.
    # The pattern: [CH2](O)[C;H]([CH2](O))O[C](=O)[*]
    #  • [CH2](O): matches a CH2 bearing an OH group (terminal glycerol arm).
    #  • [C;H]: matches the central glycerol carbon (with one hydrogen).
    #  • ([CH2](O)) : the other terminal glycerol arm.
    #  • O[C](=O)[*]: represents the ester linkage to an acyl group.
    smarts = "[CH2](O)[C;H]([CH2](O))O[C](=O)[*]"
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        return False, "Failed to parse SMARTS pattern for 2-monoglyceride"
    
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "No valid 2-monoglyceride fragment found"
    
    # For each match, check the acyl chain length (heuristic: require at least 5 contiguous aliphatic carbons)
    # In the SMARTS, mapping is as follows:
    #   idx0: terminal CH2(O)
    #   idx1: central glycerol carbon
    #   idx2: second terminal CH2(O)
    #   idx3: ester oxygen
    #   idx4: carbonyl carbon
    #   idx5: first acyl atom (start of the acyl chain)
    for match in matches:
        if len(match) < 6:
            continue  # not enough mapped atoms
        
        acyl_start_idx = match[5]
        acyl_start = mol.GetAtomWithIdx(acyl_start_idx)
        
        # Perform a breadth-first search from the acyl start atom to count contiguous aliphatic carbon atoms.
        # Exclude the carbonyl carbon (the neighbor that is part of the ester) to avoid going backwards.
        visited = set()
        queue = [acyl_start_idx]
        acyl_carbons = 0
        
        # Get the ester carbonyl atom index to exclude that branch.
        carbonyl_idx = match[4]
        
        while queue:
            current_idx = queue.pop(0)
            if current_idx in visited:
                continue
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            # Only count aliphatic carbons (atomic number 6, not aromatic)
            if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
                acyl_carbons += 1
            # Continue traversing neighbors if they are aliphatic carbons. Do not go back into the ester linkage.
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in visited:
                    continue
                if nbr_idx == carbonyl_idx:
                    continue
                # We only traverse into atoms that are part of an aliphatic chain (carbon atoms)
                if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                    queue.append(nbr_idx)
        if acyl_carbons < 5:
            # This branch's acyl chain is too short; try next match.
            continue
        
        # If a match with a sufficiently long acyl chain is found, return True.
        return True, "Contains glycerol backbone with an acyl ester at the 2-position (2-monoglyceride)"
    
    # If none of the matches qualify (e.g., acyl chain too short)
    return False, "2-monoglyceride fragment found but acyl chain is too short"

# Example usage (for testing purposes):
if __name__ == '__main__':
    test_examples = [
        "O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # MG(0:0/22:6(...)/0:0)
        "C(\\C[C@H]1[C@H](CC([C@@H]1/C=C/[C@H](CCCCC)O)=O)O)=C\\CCCC(=O)OC(CO)CO",  # prostaglandin D2 2-glyceryl ester
        "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)C(CO)CO",  # MG(0:0/24:6(...)/0:0)
        "[H]C(CO)(CO)OC(=O)CCC/C=C/C/C=C/C/C=C/C/C=C/CCCCC",  # 2-arachidonoylglycerol
        "CCO",  # obviously not a 2-monoglyceride
    ]
    
    for s in test_examples:
        valid, reason = is_2_monoglyceride(s)
        print(f"SMILES: {s}\nResult: {valid}\nReason: {reason}\n")