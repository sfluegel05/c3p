"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Straight-chain saturated fatty acid
Definition: Any saturated fatty acid lacking a side-chain.
This implementation requires that the molecule has exactly one carboxyl group and that 
all carbon atoms (ignoring heteroatoms and substituents that are not carbon) form a single, unbranched chain.
"""

from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid.
    The molecule must contain exactly one carboxyl group ([CX3](=O)[OX1H0,OX2H1]) 
    and the carbon skeleton (using only C–C single bonds) must form one linear chain (i.e. no branching).
    Non-carbon atoms (even if bound to carbons along the main chain) are ignored.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a straight-chain saturated fatty acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS pattern for a carboxyl group.
    # The pattern matches a carbon with a double-bonded oxygen and an -OH, -O– or similar.
    carboxyl_smarts = "[CX3](=O)[OX1H0,OX2H1]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    # We expect exactly one carboxyl group.
    if len(carboxyl_matches) != 1:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, expected exactly 1"
    
    # The first atom in the match is taken as the carboxyl carbon.
    carboxyl_idx = carboxyl_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    
    # Build a connectivity graph for carbon atoms (atomic num 6) via single bonds only.
    # We ignore bonds to non-carbon atoms (e.g., O, H, etc).
    carbon_indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_indices.append(atom.GetIdx())
    
    # Create graph as a dictionary: key = carbon atom index; value = list of neighboring carbon indices (via single bonds).
    carbon_graph = {idx: [] for idx in carbon_indices}
    for bond in mol.GetBonds():
        # Only consider bonds that are single and connect two carbon atoms.
        if bond.GetBondType() != Chem.BondType.SINGLE:
            # If there is any carbon–carbon bond that is not a single bond, the fatty acid is not saturated.
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                return False, "Carbon–carbon unsaturated bond found in the molecule"
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            # Add undirected edge.
            idx1, idx2 = a1.GetIdx(), a2.GetIdx()
            # Only add bond if it is single.
            if bond.GetBondType() == Chem.BondType.SINGLE:
                carbon_graph[idx1].append(idx2)
                carbon_graph[idx2].append(idx1)
    
    # Check that the carboxyl carbon is in our carbon graph.
    if carboxyl_idx not in carbon_graph:
        return False, "Carboxyl carbon not found in carbon connectivity graph"
    
    # Next, ensure that the carbon subgraph is connected.
    # Do a simple breadth-first search.
    visited = set()
    stack = [carboxyl_idx]
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        for nbr in carbon_graph[current]:
            if nbr not in visited:
                stack.append(nbr)
    total_carbons = len(carbon_graph)
    if len(visited) != total_carbons:
        return False, "Extra carbon atoms outside the main carbon chain were found (graph not connected)"
    
    # Now check that the carbon connectivity forms a linear (unbranched) chain.
    # In a straight chain (with >1 atom), exactly two carbons appear as end-points (degree 1) and the rest as interior (degree 2).
    # If only one carbon exists then degree 0 is acceptable.
    degrees = {idx: len(neighbors) for idx, neighbors in carbon_graph.items()}
    # For a fatty acid, we need at least the carboxyl carbon and one more carbon.
    if total_carbons < 2:
        return False, "Not enough carbon atoms to form a fatty acid chain"
    
    end_points = [idx for idx, deg in degrees.items() if deg == 1]
    # Also, if any carbon has degree 0 (isolated) or degree >2, then there is branching.
    for idx, deg in degrees.items():
        if deg == 0:
            return False, "Isolated carbon atom found"
        if deg > 2:
            return False, f"Carbon atom {idx} has degree {deg}, indicating branching or extra substituents"
    
    if len(end_points) != 2:
        return False, f"Expected 2 chain end-points; found {len(end_points)} (chain might be cyclic or branched)"
    
    # Finally, check that the carboxyl carbon is one of the endpoints or connected in the chain.
    # (In a fatty acid, the COOH carbon is usually terminal.)
    if degrees[carboxyl_idx] != 1:
        return False, "Carboxyl carbon is not at the chain terminus"
    
    return True, "Molecule is a straight-chain saturated fatty acid with no side-chain substituents"

# Example usage: (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
        "CC(=O)CCC(O)=O",                        # 4-oxopentanoic acid (should be True)
        "CCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",        # heptacosanoic acid (should be True)
        "CCCC(CCCCCCCCCCCCCC(O)=O)O",              # example with branching (should be False)
        "C(C[N+](C)(C)C)OP(=O)([O-])OC[C@H](OC)COC(CCCCCCCCCCCCCCC)=O",  # phospholipid; missing free COOH (False)
    ]
    for smi in test_smiles:
        valid, reason = is_straight_chain_saturated_fatty_acid(smi)
        print(f"SMILES: {smi} -> {valid}: {reason}")