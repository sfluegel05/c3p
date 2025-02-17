"""
Classifies: CHEBI:17389 2-monoglyceride
"""
#!/usr/bin/env python3
"""
Classifies: 2-monoglyceride
A 2-monoglyceride is defined as a monoglyceride in which the acyl substituent is located at position 2.
That is, the molecule should contain a glycerol backbone of the form:
   HOCH2–CH(OC(=O)R)–CH2OH
We first search for the backbone using a SMARTS pattern that requires that none of the three backbone heavy atoms are in a ring.
Then, we find the unique ester branch from the central carbon and check that its acyl substituent is aliphatic and sufficiently long.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    
    A 2-monoglyceride contains a glycerol backbone where the acyl (ester) group is attached at 
    the secondary (2-) position, i.e. the central –OH is replaced by –O–C(=O)R.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (to ensure –OH groups are explicitly seen)
    mol_with_H = Chem.AddHs(mol)
    
    # Define the SMARTS for the glycerol backbone.
    # We require that each backbone heavy carbon is not in a ring (!R).
    # The pattern captures: a primary CH2 bearing an –OH, connected to a secondary CH (with one hydrogen)
    # that has an -O substituent that in turn is attached to a carbonyl, followed by another primary CH2 bearing an –OH.
    pattern_smarts = ("[CH2;H2;!R]([OX2H])-"  # terminal CH2OH
                      "[C;H1;!R](O[C;X3](=O)[*])-"  # central CH substituted with an ester branch
                      "[CH2;H2;!R]([OX2H])")  # terminal CH2OH
    pattern = Chem.MolFromSmarts(pattern_smarts)
    if pattern is None:
        return False, "Failed to generate SMARTS pattern"
    
    backbone_matches = mol_with_H.GetSubstructMatches(pattern)
    if not backbone_matches:
        return False, "No 2-monoglyceride (sn-2 ester on glycerol) fragment found"
    
    # In a monoglyceride, expect exactly one glycerol backbone fragment.
    if len(backbone_matches) > 1:
        return False, f"Found {len(backbone_matches)} glycerol backbone fragments; expected exactly one for a monoglyceride"
    
    # The match tuple should represent the three backbone carbons.
    # match[0] -> first CH2; match[1] -> central CH; match[2] -> last CH2.
    match = backbone_matches[0]
    central_idx = match[1]
    central_atom = mol_with_H.GetAtomWithIdx(central_idx)
    
    # Identify the unique ester oxygen substituent attached to the central carbon.
    # The central carbon in glycerol normally has 3 heavy-atom neighbors:
    #   two backbone carbons (from the match) and one oxygen attached as the ester.
    backbone_set = set(match)
    ester_oxys = [nbr for nbr in central_atom.GetNeighbors() 
                  if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in backbone_set]
    if len(ester_oxys) != 1:
        return False, "Central carbon does not have a unique ester oxygen substituent"
    ester_oxy = ester_oxys[0]
    
    # From the ester oxygen, find the carbonyl carbon.
    # Typically, the oxygen is connected to a carbon bearing a double bond to oxygen.
    neighbors_oxy = [nbr for nbr in ester_oxy.GetNeighbors() if nbr.GetIdx() != central_idx and nbr.GetAtomicNum() == 6]
    if not neighbors_oxy:
        return False, "Ester oxygen does not connect to any carbonyl carbon"
    
    carbonyl = None
    for nbr in neighbors_oxy:
        bond = mol_with_H.GetBondBetweenAtoms(ester_oxy.GetIdx(), nbr.GetIdx())
        if bond is not None and bond.GetBondType() == rdchem.BondType.DOUBLE:
            carbonyl = nbr
            break
    if carbonyl is None:
        return False, "No carbonyl bond found from ester oxygen"
    
    # Identify the acyl substituent attached to the carbonyl carbon.
    # Exclude the ester oxygen; the remaining neighbor is our acyl group.
    acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors() if nbr.GetIdx() != ester_oxy.GetIdx()]
    if not acyl_neighbors:
        return False, "No acyl substituent attached to the carbonyl"
    acyl_start = acyl_neighbors[0]
    
    # Now verify that the acyl group is aliphatic (non-aromatic) and sufficiently long.
    # We perform a breadth-first traversal from the acyl substituent,
    # excluding any atoms that are part of the backbone (central ester, carbonyl, ester oxygen).
    excluded = {central_atom.GetIdx(), ester_oxy.GetIdx(), carbonyl.GetIdx()}
    queue = [acyl_start.GetIdx()]
    visited = set()
    acyl_carbon_count = 0
    while queue:
        current_idx = queue.pop(0)
        if current_idx in visited:
            continue
        visited.add(current_idx)
        atom = mol_with_H.GetAtomWithIdx(current_idx)
        # Count only carbon atoms in the acyl chain.
        if atom.GetAtomicNum() == 6:
            acyl_carbon_count += 1
            if atom.GetIsAromatic():
                return False, "Acyl substituent is aromatic, not a fatty acid chain"
        # Add neighbors that are not excluded and not already visited.
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in visited and nbr.GetIdx() not in excluded:
                queue.append(nbr.GetIdx())
    
    if acyl_carbon_count < 5:
        return False, f"Acyl substituent too short, only {acyl_carbon_count} carbon atoms found"
    
    return True, "Contains a glycerol backbone with an acyl ester at the 2-position (2-monoglyceride)"

# Example usage (for testing):
if __name__ == "__main__":
    # Testing on a few examples (true positives and some negatives)
    test_smiles = [
        "O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # MG(0:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0)
        "C(\\C[C@H]1[C@H](CC([C@@H]1/C=C/[C@H](CCCCC)O)=O)O)=C\\CCCC(=O)OC(CO)CO",  # prostaglandin D2 2-glyceryl ester
        "CCCCCCCCCCCCCCCC(=O)OC(CO)CO",  # 2-palmitoylglycerol (should be false, if the acyl chain does not qualify)
        "CCO",  # Clearly not a 2-monoglyceride
    ]
    
    for s in test_smiles:
        valid, reason = is_2_monoglyceride(s)
        print(f"SMILES: {s}\nResult: {valid}, Reason: {reason}\n")