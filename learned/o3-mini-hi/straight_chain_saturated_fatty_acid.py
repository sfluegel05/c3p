"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Straight-chain saturated fatty acid
Definition: Any saturated fatty acid lacking a side-chain.
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid.
    It checks that the molecule contains exactly one carboxyl group and that
    the remainder of the carbon skeleton is a single, unbranched, saturated (i.e. single bonds)
    chain. It ignores non-carbon atoms (e.g. hydroxyl groups) in the chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a straight-chain saturated fatty acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern to capture a carboxyl group.
    # This pattern matches a carbon with a double bonded oxygen and an –OH or –O– (possibly deprotonated).
    carboxyl_smarts = "[CX3](=O)[OX1H0-,OX2H1]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    # We expect exactly one carboxyl group.
    if len(carboxyl_matches) != 1:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, expected exactly 1"
    
    # The first atom in the match is the carbonyl carbon.
    carboxyl_idx = carboxyl_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    
    # Identify the α-carbon (the carbon directly bonded to the carboxyl carbon).
    alpha_carbon = None
    for neighbor in carboxyl_atom.GetNeighbors():
        # The carboxyl carbon should be bonded to oxygen(s) and exactly one carbon.
        if neighbor.GetAtomicNum() == 6:
            alpha_carbon = neighbor
            break
    if alpha_carbon is None:
        return False, "No α-carbon found adjacent to the carboxyl group"

    # We now traverse the carbon chain starting from the α-carbon.
    # We must ensure that (a) all bonds along the chain are single, and (b) no branching (extra carbon substituents) occurs.
    chain_atoms = set()  # store indices of carbon atoms in the chain (excluding the carboxyl carbon)
    visited = set()
    branching_found = False

    def traverse(atom, parent):
        nonlocal branching_found
        idx = atom.GetIdx()
        if idx in visited:
            return
        visited.add(idx)
        chain_atoms.add(idx)
        # Get carbon neighbors (ignore the parent).
        carbon_neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == parent.GetIdx():
                continue
            if nbr.GetAtomicNum() == 6:
                carbon_neighbors.append(nbr)
        # If there is more than one carbon neighbor, then there is a branch.
        if len(carbon_neighbors) > 1:
            branching_found = True
            return
        # Traverse further along the chain.
        for nbr in carbon_neighbors:
            bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
            # Check that the C–C bond is single.
            if bond.GetBondType() != Chem.BondType.SINGLE:
                branching_found = True
                return
            traverse(nbr, atom)
    
    traverse(alpha_carbon, carboxyl_atom)
    if branching_found:
        return False, "Chain is branched or contains unsaturated carbon–carbon bonds"
    
    # Now verify that the only carbon atoms in the molecule belong to the main chain.
    all_carbon_indices = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6}
    # Expected carbons: the carboxyl carbon plus all carbons encountered in the chain traversal.
    expected_carbons = set(chain_atoms) | {carboxyl_idx}
    if all_carbon_indices != expected_carbons:
        return False, "Extra carbon atoms outside the main chain were found (possible side-chain substituents)"
    
    return True, "Molecule is a straight-chain saturated fatty acid with no side-chain substituents"

# Example usage (you can remove or comment out these examples when integrating into your system):
if __name__ == "__main__":
    test_smiles = [
        "CCCC(O)=O",               # butyric acid (should be True)
        "CC(=O)CCC(O)=O",          # 4-oxopentanoic acid (containing a carbonyl internally – not straight saturated chain)
        "CCCC(C)C(O)=O"            # example with branching (should be False)
    ]
    for smi in test_smiles:
        valid, reason = is_straight_chain_saturated_fatty_acid(smi)
        print(f"SMILES: {smi} -> {valid}: {reason}")