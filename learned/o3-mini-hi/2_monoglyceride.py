"""
Classifies: CHEBI:17389 2-monoglyceride
"""
#!/usr/bin/env python3
"""
Classifies: 2-monoglyceride
A 2-monoglyceride is defined as a monoglyceride in which the acyl substituent is located at the 2-position.
That is, the glycerol backbone is modified from HOCH2–CH(OH)–CH2OH into HOCH2–CH(O–C(=O)R)–CH2OH.
This program uses a heuristic approach to first identify a candidate “central” carbon (sn-2)
by looking for a sp3 carbon connected to exactly three heavy neighbors: two carbons (which should lead to CH2OH groups)
and one oxygen (whose bond should be part of an ester group, i.e. connected to a carbonyl).
Then it verifies that the acyl chain attached at that ester bond is aliphatic and of a sufficient length.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    
    A 2-monoglyceride contains a glycerol backbone where the acyl (ester) group is attached at 
    the central (2-) position. That is, one of the –OH groups of glycerol is replaced by an acyl ester.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to reliably see OH groups, etc.
    mol = Chem.AddHs(mol)
    
    # Loop over all atoms in the molecule to identify a candidate “central” carbon.
    # In a 2-monoglyceride, the central carbon should be sp3, not in a ring,
    # and it should be connected to exactly three heavy (non-hydrogen) atoms:
    # two carbons (which should lead to CH2OH groups) and one oxygen which must be part of an ester.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # consider only carbon atoms
        
        # We require that the candidate carbon is not in a ring.
        if atom.IsInRing():
            continue
        
        # Count heavy-atom neighbors (exclude hydrogens).
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 3:
            continue
        
        # Among the heavy neighbors, exactly one must be oxygen (ester oxygen) and two carbons.
        oxy_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetAtomicNum() == 8]
        carbon_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetAtomicNum() == 6]
        if len(oxy_neighbors) != 1 or len(carbon_neighbors) != 2:
            continue
        
        ester_oxygen = oxy_neighbors[0]
        
        # Check that the two remaining neighbors (the “terminal” groups) are likely CH2OH groups.
        # We will look that each such carbon has at least one oxygen neighbor (an –OH) that is not our central carbon.
        terminals_ok = True
        for term in carbon_neighbors:
            term_oxygens = [nbr for nbr in term.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != atom.GetIdx()]
            if not term_oxygens:
                terminals_ok = False
                break
        if not terminals_ok:
            continue
        
        # Now check the connectivity of the ester oxygen.
        # It should be connected to our central carbon and also to a carbonyl carbon (i.e. a carbon with a double bond to an oxygen).
        # Get the neighbor(s) of the ester oxygen other than the candidate central carbon.
        neighboring_carbons = [nbr for nbr in ester_oxygen.GetNeighbors() if nbr.GetIdx() != atom.GetIdx() and nbr.GetAtomicNum() == 6]
        if not neighboring_carbons:
            continue
        carbonyl = None
        for nbr in neighboring_carbons:
            bond = mol.GetBondBetweenAtoms(ester_oxygen.GetIdx(), nbr.GetIdx())
            if bond is not None and bond.GetBondType() == rdchem.BondType.DOUBLE:
                carbonyl = nbr
                break
        if carbonyl is None:
            continue  # ester oxygen is not bonded to a carbonyl carbon
        
        # From the carbonyl carbon, identify the acyl substituent.
        acyl_neighbors = [nbr for nbr in carbonyl.GetNeighbors() if nbr.GetIdx() != ester_oxygen.GetIdx()]
        if not acyl_neighbors:
            continue
        acyl_start = acyl_neighbors[0]
        
        # Traverse the acyl chain using a breadth-first search.
        # We count the number of carbon atoms in the acyl chain. We also check if any atom is aromatic,
        # which would rule it out as a "fatty acid" chain.
        excluded = {atom.GetIdx(), ester_oxygen.GetIdx(), carbonyl.GetIdx()}
        queue = [acyl_start.GetIdx()]
        visited = set()
        acyl_carbon_count = 0
        while queue:
            cur_idx = queue.pop(0)
            if cur_idx in visited:
                continue
            visited.add(cur_idx)
            current_atom = mol.GetAtomWithIdx(cur_idx)
            if current_atom.GetAtomicNum() == 6:
                acyl_carbon_count += 1
                if current_atom.GetIsAromatic():
                    return False, "Acyl substituent is aromatic, not an aliphatic fatty acid chain"
            # Add neighbors that are not in the excluded set.
            for nbr in current_atom.GetNeighbors():
                if nbr.GetIdx() not in visited and nbr.GetIdx() not in excluded:
                    queue.append(nbr.GetIdx())
        
        # Heuristically, require the acyl chain to have at least 5 carbon atoms.
        if acyl_carbon_count < 5:
            # Too short to be a typical fatty acyl chain.
            continue
        
        # If we reach here, we have found a candidate: a central glycerol carbon (with two terminal CH2OH groups)
        # having a unique oxygen substituent leading via a carbonyl to an aliphatic acyl chain.
        return True, "Contains glycerol backbone with an acyl ester at the 2-position (2-monoglyceride)"
    
    # If no candidate fragment fits, return False.
    return False, "No valid 2-monoglyceride fragment found"

# Example usage (for testing purposes):
if __name__ == '__main__':
    test_examples = [
        # These are some of the provided examples:
        "O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # MG(0:0/22:6(...)/0:0)
        "C(\\C[C@H]1[C@H](CC([C@@H]1/C=C/[C@H](CCCCC)O)=O)O)=C\\CCCC(=O)OC(CO)CO",  # prostaglandin D2 2-glyceryl ester
        "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)C(CO)CO",  # MG(0:0/24:6(...)/0:0)
        "[H]C(CO)(CO)OC(=O)CCC/C=C/C/C=C/C/C=C/C/C=C/CCCCC",  # 2-arachidonoylglycerol
        # An obviously non-2-monoglyceride example:
        "CCO",  
    ]
    
    for s in test_examples:
        valid, reason = is_2_monoglyceride(s)
        print(f"SMILES: {s}\nResult: {valid}\nReason: {reason}\n")