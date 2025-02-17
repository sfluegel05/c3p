"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: Monoacylglycerol
Definition: A glyceride in which exactly one of the three –OH groups of a glycerol backbone is esterified
while the other two remain as free hydroxyl groups.
Heuristic approach:
  1. Parse the SMILES string.
  2. Count the total number of ester groups using a generic ester pattern "[OX2][CX3](=O)"
     – a monoacylglycerol must have exactly 1 such ester.
  3. Count free hydroxyl groups using "[OX2H]".
     – we require exactly 2 free –OH groups.
  4. Identify the candidate glycerol backbone carbons:
       • For each free –OH, take its attached carbon.
       • For the ester group, take the oxygen that is not the carbonyl and get its attached carbon.
     We expect these to produce exactly 3 unique carbon atoms.
  5. Check that these three candidate carbons are “contiguous” in the molecule.
     In a glycerol backbone they must form a linear chain: two of them are bonded directly to the central one
     (thus giving two bonds overall, with degrees d=1, d=2, d=1).
  6. If all are met, we classify the molecule as a monoacylglycerol.
  
If any check fails, a reason is reported.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string using refined heuristics.
    
    Approach:
     - Ensure exactly one ester (acyl) group is present.
     - Ensure exactly two free hydroxyl groups are present.
     - Identify three candidate glycerol carbons:
         • Two coming from free –OH groups
         • One coming from the alcohol oxygen that is esterified.
     - Check that these three candidate carbons form a contiguous (linear, acyclic) 3-carbon backbone.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a monoacylglycerol, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Check for exactly one ester group.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Number of ester groups is {len(ester_matches)}; expected exactly 1 for a monoacylglycerol"
    
    # Step 2. Check for exactly two free hydroxyl groups.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if len(oh_matches) != 2:
        return False, f"Number of free hydroxyl groups is {len(oh_matches)}; expected exactly 2 for a monoacylglycerol"
    
    # Step 3. Identify candidate glycerol carbons.
    # From the free hydroxyl groups, take the attached carbon (if any)
    candidate_carbons = set()
    for match in oh_matches:
        # match is a tuple with one oxygen index
        oh_atom = mol.GetAtomWithIdx(match[0])
        # Get neighbor carbons
        carbon_neighbors = [nbr.GetIdx() for nbr in oh_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not carbon_neighbors:
            continue
        # (In a typical -OH, there should be one bonded carbon.)
        candidate_carbons.update(carbon_neighbors)
    
    # From the ester group match, determine the oxygen that connects glycerol to the acyl chain.
    # The match tuple is (o_idx, c_idx) where the O is the alcohol oxygen (should be on glycerol)
    ester_match = ester_matches[0]
    ester_o = mol.GetAtomWithIdx(ester_match[0])
    ester_candidate = None
    # The ester oxygen is typically connected to the carbonyl carbon and to the glycerol carbon.
    for nbr in ester_o.GetNeighbors():
        # Exclude the carbonyl carbon (given by ester_match[1]); the other carbon is the candidate.
        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != ester_match[1]:
            ester_candidate = nbr.GetIdx()
            break
    if ester_candidate is not None:
        candidate_carbons.add(ester_candidate)
    
    # We expect exactly three distinct candidate carbons.
    candidate_carbons = list(candidate_carbons)
    if len(candidate_carbons) != 3:
        return False, f"Found {len(candidate_carbons)} candidate glycerol carbons; expected 3 (from 2 free OH and 1 ester linkage)"
    
    # Step 4. Check connectivity of the candidate carbons.
    # In a glycerol backbone, the three carbons are connected in a linear chain.
    # We build a simple graph: each candidate carbon is a node and we add an edge if they are directly bonded.
    edges = {idx: [] for idx in candidate_carbons}
    for i in range(len(candidate_carbons)):
        for j in range(i+1, len(candidate_carbons)):
            # Check if there is a direct bond between candidate carbons.
            bond = mol.GetBondBetweenAtoms(candidate_carbons[i], candidate_carbons[j])
            if bond is not None:
                edges[candidate_carbons[i]].append(candidate_carbons[j])
                edges[candidate_carbons[j]].append(candidate_carbons[i])
    
    # For a linear three-atom chain the edge counts should be: one carbon with 2 bonds (middle) and two carbons with 1 each.
    total_edges = sum(len(neighs) for neighs in edges.values())
    if total_edges != 4:
        return False, "Candidate glycerol carbons are not properly connected (expected 2 bonds total)"
    
    degrees = [len(edges[idx]) for idx in candidate_carbons]
    if sorted(degrees) != [1, 1, 2]:
        return False, "Candidate glycerol carbons do not form a linear (contiguous) 3-carbon chain"
    
    # (Optional) Check that molecular weight is above a reasonable threshold.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 140:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a monoacylglycerol"
    
    return True, "Contains a contiguous 3‐carbon (glycerol) backbone with exactly one acyl (ester) substitution and two free hydroxyl groups, consistent with a monoacylglycerol"


# For testing (this section will run if the file is executed directly)
if __name__ == "__main__":
    # Examples to test (feel free to add more)
    test_smiles_list = [
        "O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # MG(0:0/22:6... example
        "CCCCCCCC(=O)OC[C@@H](O)CO",  # 1-octanoyl-sn-glycerol
        "O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)C(CO)CO",  # MG(0:0/24:6... example
    ]
    for s in test_smiles_list:
        result, reason = is_monoacylglycerol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")